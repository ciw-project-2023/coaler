//
// Created by malte on 12/4/23.
//

#include "Matcher.hpp"

#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <spdlog/spdlog.h>

#include <queue>

#include "GraphMol/ChemTransforms/ChemTransforms.h"
#include "GraphMol/DistGeomHelpers/Embedder.h"
#include "GraphMol/FMCS/FMCS.h"

namespace coaler::core {

    void Matcher::murckoPruningRecursive(RDKit::RWMOL_SPTR mol, int atomID, int parentID, std::vector<bool>& visit,
                                         std::vector<int>& delAtoms, std::vector<std::pair<int, int>>& delBonds) {
        // mark atom as visited and find all neighbor atoms of atom with atomID
        visit.at(atomID) = true;
        for (const auto& nbri : boost::make_iterator_range(mol->getAtomNeighbors(mol->getAtomWithIdx(atomID)))) {
            // only visit atoms not yet visited
            if (nbri == parentID || visit.at(nbri)) {
                continue;
            }
            murckoPruningRecursive(mol, nbri, atomID, visit, delAtoms, delBonds);

            // find neighbors that are not on the to-be-deleted-list and count them
            int numNbrs = 0;
            for (const auto& nbriNext : boost::make_iterator_range(mol->getAtomNeighbors(mol->getAtomWithIdx(nbri)))) {
                const auto& nbrNext = (*mol)[nbriNext];
                if (std::find(delAtoms.begin(), delAtoms.end(), nbriNext) == delAtoms.end()) {
                    numNbrs++;
                }
            }

            // atoms with only one neighbor (parent) will be deleted after recursive call.
            // They are not deleted here cause internal IDs are updated after deletion.
            if (numNbrs == 1) {
                delBonds.push_back(std::make_pair(nbri, atomID));
                delAtoms.push_back(nbri);
            }
        }
    }

    std::optional<CoreResult> Matcher::calculateCoreMcs(RDKit::MOL_SPTR_VECT mols, int numOfThreads) {
        // Generates all parameters needed for RDKit::findMCS()
        RDKit::MCSParameters mcsParams;
        RDKit::MCSAtomCompareParameters atomCompParams;
        atomCompParams.MatchChiralTag = true;
        atomCompParams.MatchFormalCharge = true;
        atomCompParams.MatchIsotope = false;
        atomCompParams.MatchValences = true;
        atomCompParams.RingMatchesRingOnly = true;
        atomCompParams.CompleteRingsOnly = false;
        mcsParams.AtomCompareParameters = atomCompParams;

        RDKit::MCSBondCompareParameters bondCompParams;
        bondCompParams.MatchStereo = true;
        bondCompParams.RingMatchesRingOnly = false;
        bondCompParams.CompleteRingsOnly = true;
        bondCompParams.MatchFusedRings = true;
        bondCompParams.MatchFusedRingsStrict = false;
        mcsParams.BondCompareParameters = bondCompParams;

        mcsParams.setMCSAtomTyperFromEnum(RDKit::AtomCompareAnyHeavyAtom);
        mcsParams.setMCSBondTyperFromEnum(RDKit::BondCompareAny);

        RDKit::MCSResult const mcs = RDKit::findMCS(mols, &mcsParams);
        if (mcs.QueryMol == nullptr) {
            return std::nullopt;
        }

        spdlog::info("MCS: {}", mcs.SmartsString);

        RDKit::RWMol first = *mols.at(0);

        RDKit::DGeomHelpers::EmbedParameters params;
        RDKit::DGeomHelpers::EmbedMolecule(first, params);

        std::vector<std::pair<int, double>> result;
        RDKit::MMFF::MMFFOptimizeMoleculeConfs(first, result);

        RDKit::SubstructMatchParameters substructMatchParams;
        substructMatchParams.useChirality = true;
        substructMatchParams.useEnhancedStereo = true;
        substructMatchParams.aromaticMatchesConjugated = true;
        substructMatchParams.numThreads = numOfThreads;

        std::vector<RDKit::MatchVectType> structMatches
            = RDKit::SubstructMatch(first, *mcs.QueryMol, substructMatchParams);
        assert(!structMatches.empty());

        AtomMap moleculeCoreCoords;
        RDKit::Conformer conformer = first.getConformer(0);
        for (const auto& [queryId, molId] : structMatches.at(0)) {
            const RDGeom::Point3D atomCoords = conformer.getAtomPos(molId);
            moleculeCoreCoords.emplace(queryId, atomCoords);
        }

        return std::make_pair(mcs.QueryMol, moleculeCoreCoords);
    }

    std::optional<CoreResult> Matcher::calculateCoreMurcko(RDKit::MOL_SPTR_VECT mols, int numOfThreads) {
        auto mcs = Matcher::calculateCoreMcs(mols, numOfThreads);
        if (!mcs.has_value()) {
            return std::nullopt;
        }

        RDKit::RWMol mcsRWMol = *mcs.value().first;
        RDKit::MolOps::sanitizeMol(mcsRWMol);
        if (mcsRWMol.getRingInfo()->numRings() == 0) {
            return std::nullopt;
        }

        std::vector<bool> visit(mcsRWMol.getNumAtoms(), false);
        std::vector<int> atomQueue;

        std::vector<std::vector<int>> ringVec = mcsRWMol.getRingInfo()->atomRings();
        for (auto ring : ringVec) {
            for (int atomID : ring) {
                if (visit.at(atomID)) {
                    continue;
                }
                atomQueue.push_back(atomID);
            }
        }

        std::vector<int> delAtoms;
        std::vector<std::pair<int, int>> delBonds;
        RDKit::RWMOL_SPTR mcsRWPtr = boost::make_shared<RDKit::RWMol>(mcsRWMol);
        for (auto atom : atomQueue) {
            for (const auto& nbri :
                 boost::make_iterator_range(mcsRWMol.getAtomNeighbors(mcsRWMol.getAtomWithIdx(atom)))) {
                const auto& nbr = (mcsRWMol)[nbri];
                murckoPruningRecursive(mcsRWPtr, nbri, atom, visit, delAtoms, delBonds);
            }
        }

        std::sort(delAtoms.begin(), delAtoms.end(), std::greater<int>());
        for (auto [atom1, atom2] : delBonds) {
            mcsRWPtr->removeBond(atom1, atom2);
        }
        for (auto atom : delAtoms) {
            mcsRWPtr->removeAtom(atom);
        }

        spdlog::info("Murco: {}", RDKit::MolToSmarts(*mcsRWPtr));

        RDKit::RWMol first = *mols.at(0);

        RDKit::DGeomHelpers::EmbedParameters params;
        RDKit::DGeomHelpers::EmbedMolecule(first, params);

        std::vector<std::pair<int, double>> result;
        RDKit::MMFF::MMFFOptimizeMoleculeConfs(first, result);

        RDKit::SubstructMatchParameters substructMatchParams;
        substructMatchParams.useChirality = true;
        substructMatchParams.useEnhancedStereo = true;
        substructMatchParams.aromaticMatchesConjugated = true;
        substructMatchParams.numThreads = numOfThreads;

        std::vector<RDKit::MatchVectType> structMatches = RDKit::SubstructMatch(first, *mcsRWPtr, substructMatchParams);
        assert(!structMatches.empty());

        AtomMap moleculeCoreCoords;
        RDKit::Conformer conformer = first.getConformer(0);
        for (const auto& [queryId, molId] : structMatches.at(0)) {
            const RDGeom::Point3D atomCoords = conformer.getAtomPos(molId);
            moleculeCoreCoords.emplace(queryId, atomCoords);
        }

        return std::make_pair(mcsRWPtr, moleculeCoreCoords);
    }

}  // namespace coaler::core