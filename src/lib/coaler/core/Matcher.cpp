//
// Created by malte on 12/4/23.
//

#include "Matcher.hpp"

#include <omp.h>
#include <spdlog/spdlog.h>

#include "GraphMol/ChemTransforms/ChemTransforms.h"
#include "GraphMol/DistGeomHelpers/Embedder.h"
#include "GraphMol/FMCS/FMCS.h"
#include "GraphMol/ForceFieldHelpers/MMFF/MMFF.h"
#include "GraphMol/RWMol.h"
#include "GraphMol/SmilesParse/SmartsWrite.h"
#include "GraphMol/SmilesParse/SmilesWrite.h"

namespace coaler::core {
    Matcher::Matcher(int threads) : m_threads(threads) {}

    std::optional<CoreResult> Matcher::calculateCoreMcs(RDKit::MOL_SPTR_VECT &mols) {
        // Generates all parameters needed for RDKit::findMCS()
        auto mcsParams = Matcher::getRelaxedMCSParams();
        RDKit::MCSResult const globalMcs = RDKit::findMCS(mols, &mcsParams);
        if (globalMcs.QueryMol == nullptr) {
            return std::nullopt;
        }

        spdlog::info("MCS: {}", globalMcs.SmartsString);

        auto pairwiseMcs = this->getPairwiseMCS(mols);
        auto [refId, reference] = this->findReferenceMolecule(mols, pairwiseMcs);
        auto ref = this->buildMolConformerForQuery(reference, *globalMcs.QueryMol);

        auto matches = RDKit::SubstructMatch(*ref, *globalMcs.QueryMol, this->getMatchParams());
        assert(!matches.empty());

        auto match = matches.at(0);
        std::unordered_map<int, int> coreToRef;
        for (auto const &[queryId, molId] : match) {
            coreToRef[queryId] = molId;
        }

        return CoreResult{globalMcs.QueryMol, refId, ref, coreToRef, pairwiseMcs};
    }

    std::optional<CoreResult> Matcher::calculateCoreMurcko(RDKit::MOL_SPTR_VECT &mols) {
        // calculate MCS first and sanitize molecule
        auto mcs = Matcher::calculateCoreMcs(mols);
        if (!mcs.has_value()) {
            return std::nullopt;
        }

        RDKit::RWMol mcsRWMol = *mcs.value().core;
        RDKit::MolOps::sanitizeMol(mcsRWMol);
        if (mcsRWMol.getRingInfo()->numRings() == 0) {
            return std::nullopt;
        }

        // initialize vars for DFS and find all atoms that are part of any ring
        std::vector<bool> visit(mcsRWMol.getNumAtoms(), false);

        std::vector<int> ringAtoms;
        std::vector<std::vector<int>> ringVec = mcsRWMol.getRingInfo()->atomRings();
        for (const auto ring : ringVec) {
            for (const auto atomID : ring) {
                if (visit.at(atomID)) {
                    continue;
                }

                ringAtoms.push_back(atomID);
            }
        }

        // start pruning of mcs, save to be deleted atoms and bonds
        std::vector<int> delAtomsMaybe;
        std::vector<std::pair<int, int>> delBonds;
        RDKit::RWMOL_SPTR murckoPtr = boost::make_shared<RDKit::RWMol>(mcsRWMol);

        // start from each ring atom
        for (const auto atomID : ringAtoms) {
            for (int i = 0; i < murckoPtr->getNumAtoms(); i++) {
                visit.at(i) = false;
            }
            murckoPruningRecursive(murckoPtr, atomID, -1, visit, delAtomsMaybe, delBonds, ringAtoms);
        }

        // delete duplicates from deletetion lists
        sort(delAtomsMaybe.begin(), delAtomsMaybe.end());
        delAtomsMaybe.erase(unique(delAtomsMaybe.begin(), delAtomsMaybe.end()), delAtomsMaybe.end());
        sort(delBonds.begin(), delBonds.end());
        delBonds.erase(unique(delBonds.begin(), delBonds.end()), delBonds.end());

        // check atoms in delAtomMaybe if they are part of the substructure between rings. If not
        // they are added to a delAtomsDefinitely.
        // An atom is part of the substructure between at least two rings if the DFS finds two or more
        // ring atoms.
        std::vector<int> foundRingAtoms;
        std::vector<int> delAtomsDefinitely;
        for (const auto delAtomsID : delAtomsMaybe) {
            for (int i = 0; i < murckoPtr->getNumAtoms(); i++) {
                visit.at(i) = false;
            }
            foundRingAtoms.clear();

            murckoCheckDelAtoms(murckoPtr, delAtomsID, -1, visit, ringAtoms, foundRingAtoms);
            if (foundRingAtoms.size() < 2) {
                delAtomsDefinitely.push_back(delAtomsID);
            }
        }

        // remove all bonds and atoms that are not part of  the murcko scaffold.
        // Deletion of atoms needs to be in order of atomIdx (high to low) to avoid
        // deletion errors.
        std::sort(delAtomsDefinitely.begin(), delAtomsDefinitely.end(), std::greater<>());
        for (auto [atom1, atom2] : delBonds) {
            if (std::find(delAtomsDefinitely.begin(), delAtomsDefinitely.end(), atom1) != delAtomsDefinitely.end()
                && std::find(delAtomsDefinitely.begin(), delAtomsDefinitely.end(), atom2) != delAtomsDefinitely.end()) {
                murckoPtr->removeBond(atom1, atom2);
            }
        }
        for (auto atom : delAtomsDefinitely) {
            murckoPtr->removeAtom(atom);
        }

        spdlog::info("Murco: {}", RDKit::MolToSmarts(*murckoPtr));

        // Embedding of core and calculation of atomCoords
        RDKit::RWMol first = *mols.at(0);
        auto ref = this->buildMolConformerForQuery(first, *murckoPtr);

        auto matches = RDKit::SubstructMatch(*ref, *murckoPtr, this->getMatchParams());
        auto match = matches.at(0);
        std::unordered_map<int, int> coreToRef;
        for (auto const &[queryId, molId] : match) {
            coreToRef[queryId] = molId;
        }

        return CoreResult{murckoPtr, 0, ref, coreToRef};
    }

    // TODO this only creates a conformer for one specific reference point. A query can contain wildcards
    // for atoms (i.e. [#6,#7] = C or N and this can have an impact on conformers (swapping Atoms can lead
    // differnet geometries). We should take that into account and try to create a diverse set of conformers
    // thta models the output of the query more closely. For now, this is a pretty good reference though.
    RDKit::ROMOL_SPTR Matcher::buildMolConformerForQuery(RDKit::RWMol first, RDKit::ROMol query) {
        auto params = RDKit::DGeomHelpers::srETKDGv3;
        params.numThreads = m_threads;
        params.randomSeed = 42;
        params.useRandomCoords = true;
        RDKit::DGeomHelpers::EmbedMolecule(first, params);

        std::vector<std::pair<int, double>> result;
        RDKit::MMFF::MMFFOptimizeMoleculeConfs(first, result, m_threads);

        return boost::make_shared<RDKit::ROMol>(first);
    }

    void Matcher::murckoPruningRecursive(RDKit::RWMOL_SPTR &mol, int atomID, int parentID, std::vector<bool> &visit,
                                         std::vector<int> &delAtoms, std::vector<std::pair<int, int>> &delBonds,
                                         std::vector<int> &ringAtoms) {
        // mark atom as visited and find all neighbor atoms of atom with atomID
        visit.at(atomID) = true;
        for (const auto &neighborID : boost::make_iterator_range(mol->getAtomNeighbors(mol->getAtomWithIdx(atomID)))) {
            // only visit atoms not the parent, not yet visited and which are not part of any ring
            if (neighborID == parentID || visit.at(neighborID)
                || std::find(ringAtoms.begin(), ringAtoms.end(), neighborID) != ringAtoms.end()) {
                continue;
            }

            murckoPruningRecursive(mol, neighborID, atomID, visit, delAtoms, delBonds, ringAtoms);

            // find neighbors that are not on the to-be-deleted-list and count them
            int numNbrs = 0;
            for (const auto &nextNeighborID :
                 boost::make_iterator_range(mol->getAtomNeighbors(mol->getAtomWithIdx(neighborID)))) {
                if (std::find(delAtoms.begin(), delAtoms.end(), nextNeighborID) == delAtoms.end()) {
                    numNbrs++;
                }
            }

            // atoms with only one neighbor (parent) will be deleted after recursive call.
            // They are not deleted here cause internal IDs are updated after deletion.
            if (numNbrs == 1) {
                delBonds.push_back(std::make_pair(neighborID, atomID));
                delAtoms.push_back(neighborID);
            }
        }
    }

    void Matcher::murckoCheckDelAtoms(RDKit::RWMOL_SPTR &mol, int atomID, int parentID, std::vector<bool> &visit,
                                      std::vector<int> &ringAtoms, std::vector<int> &foundRingAtoms) {
        // mark atom as visited and check if they are part of any ring. If yes save and return to
        // recursive call before
        visit.at(atomID) = true;
        if (std::find(ringAtoms.begin(), ringAtoms.end(), atomID) != ringAtoms.end()) {
            foundRingAtoms.push_back(atomID);
            return;
        }
        // visit all neighbors if they are not the parent and are not yet visited
        for (const auto &neighborID : boost::make_iterator_range(mol->getAtomNeighbors(mol->getAtomWithIdx(atomID)))) {
            if (neighborID == parentID || visit.at(neighborID)) {
                continue;
            }
            murckoCheckDelAtoms(mol, neighborID, atomID, visit, ringAtoms, foundRingAtoms);
        }
    }

    /*----------------------------------------------------------------------------------------------------------------*/
    PairwiseMCSResult Matcher::getPairwiseMCS(RDKit::MOL_SPTR_VECT &mols) {
        PairwiseMCSMap strictMcsMap, relaxedMcsMap;

        omp_lock_t mapLock;
        omp_init_lock(&mapLock);

#pragma omp parallel for shared(mols, strictMcsMap, relaxedMcsMap, mapLock) default(none)
        for (auto firstId = 0; firstId < mols.size(); firstId++) {
            auto firstMol = mols.at(firstId);
            for (auto secondId = firstId + 1; secondId < mols.size(); secondId++) {
                const auto secondMol = mols.at(secondId);
                const auto key = multialign::LigandPair(firstId, secondId);

                const RDKit::MOL_SPTR_VECT molPair{firstMol, secondMol};

                const auto strictParams = getStrictMCSParams();
                const auto relaxedParams = getRelaxedMCSParams();

                const auto firstMolSmiles = RDKit::MolToSmiles(*firstMol);
                const auto secondMolSmiles = RDKit::MolToSmiles(*secondMol);

                const auto strictMcs = RDKit::findMCS(molPair, &strictParams);
                const auto relaxedMcs = RDKit::findMCS(molPair, &relaxedParams);

                auto empty = std::tuple<RDKit::MatchVectType, RDKit::MatchVectType, std::string>();
                if (strictMcs.QueryMol == nullptr) {
                    spdlog::error("no strict mcs found between {} and {}", firstMolSmiles, secondMolSmiles);
                    omp_set_lock(&mapLock);
                    strictMcsMap.emplace(key, empty);
                    omp_unset_lock(&mapLock);
                }

                if (relaxedMcs.QueryMol == nullptr) {
                    spdlog::error("no relaxed mcs found between {} and {}", firstMolSmiles, secondMolSmiles);
                    omp_set_lock(&mapLock);
                    relaxedMcsMap.emplace(key, empty);
                    omp_unset_lock(&mapLock);
                }

                auto firstStrictMatches = RDKit::SubstructMatch(*firstMol, *strictMcs.QueryMol, getMatchParams());
                auto secondStrictMatches = RDKit::SubstructMatch(*secondMol, *strictMcs.QueryMol, getMatchParams());
                auto firstRelaxedMatches = RDKit::SubstructMatch(*firstMol, *relaxedMcs.QueryMol, getMatchParams());
                auto secondRelaxedMatches = RDKit::SubstructMatch(*secondMol, *relaxedMcs.QueryMol, getMatchParams());

                if (firstStrictMatches.empty() || secondStrictMatches.empty()) {
                    spdlog::error("no strict mcs match between {} and {}", firstMolSmiles, secondMolSmiles);
                    omp_set_lock(&mapLock);
                    strictMcsMap.emplace(key, empty);
                    omp_unset_lock(&mapLock);
                }

                if (firstRelaxedMatches.empty() || secondRelaxedMatches.empty()) {
                    spdlog::error("no relaxed mcs match between {} and {}", firstMolSmiles, secondMolSmiles);
                    omp_set_lock(&mapLock);
                    relaxedMcsMap.emplace(key, empty);
                    omp_unset_lock(&mapLock);
                }

                omp_set_lock(&mapLock);
                strictMcsMap.emplace(
                    key, std::make_tuple(firstStrictMatches.at(0), secondStrictMatches.at(0), strictMcs.SmartsString));
                relaxedMcsMap.emplace(key, std::make_tuple(firstRelaxedMatches.at(0), secondRelaxedMatches.at(0),
                                                           relaxedMcs.SmartsString));
                omp_unset_lock(&mapLock);
            }
        }

        return {strictMcsMap, relaxedMcsMap};
    }

    std::pair<multialign::LigandID, RDKit::ROMol> Matcher::findReferenceMolecule(RDKit::MOL_SPTR_VECT &mols, PairwiseMCSResult &pairwiseMcsResult) {
        auto relaxedMcsMap = pairwiseMcsResult.relaxedMcsMap;

        std::map<multialign::LigandID, int> result;
        for (auto [key, value] : relaxedMcsMap) {
            auto [firstId, secondId] = key;
            auto [firstMatch, secondMatch, _] = value;

            result[firstId] += firstMatch.size();
            result[secondId] += secondMatch.size();
        }

        multialign::LigandID candidate;
        int max = 0;
        for (auto [key, value] : result) {
            if (value > max) {
                candidate = key;
                max = value;
            }
        }

        spdlog::info("reference molecule: {}", RDKit::MolToSmiles(*mols.at(candidate)));

        return {candidate, *mols.at(candidate)};
    }

    RDKit::SubstructMatchParameters Matcher::getMatchParams() const {
        RDKit::SubstructMatchParameters substructMatchParams;
        substructMatchParams.useChirality = true;
        substructMatchParams.useEnhancedStereo = true;
        substructMatchParams.aromaticMatchesConjugated = true;
        substructMatchParams.numThreads = m_threads;

        return substructMatchParams;
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    RDKit::MCSParameters Matcher::getRelaxedMCSParams() {
        RDKit::MCSParameters mcsParams;
        RDKit::MCSAtomCompareParameters atomCompParams;
        atomCompParams.MatchChiralTag = false;
        atomCompParams.MatchFormalCharge = false;
        atomCompParams.MatchIsotope = false;
        atomCompParams.MatchValences = true;
        atomCompParams.RingMatchesRingOnly = true;
        atomCompParams.CompleteRingsOnly = true;
        mcsParams.AtomCompareParameters = atomCompParams;

        RDKit::MCSBondCompareParameters bondCompParams;
        bondCompParams.MatchStereo = false;
        bondCompParams.RingMatchesRingOnly = false;
        bondCompParams.CompleteRingsOnly = true;
        bondCompParams.MatchFusedRings = true;
        bondCompParams.MatchFusedRingsStrict = false;
        mcsParams.BondCompareParameters = bondCompParams;
        mcsParams.Timeout = 1;
        mcsParams.setMCSAtomTyperFromEnum(RDKit::AtomCompareAny);
        mcsParams.setMCSBondTyperFromEnum(RDKit::BondCompareAny);
        return mcsParams;
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    RDKit::MCSParameters Matcher::getStrictMCSParams() {
        RDKit::MCSParameters mcsParams;
        RDKit::MCSAtomCompareParameters atomCompParams;
        atomCompParams.MatchChiralTag = true;
        atomCompParams.MatchFormalCharge = true;
        atomCompParams.MatchIsotope = false;
        atomCompParams.MatchValences = true;
        atomCompParams.RingMatchesRingOnly = true;
        atomCompParams.CompleteRingsOnly = true;
        mcsParams.AtomCompareParameters = atomCompParams;

        RDKit::MCSBondCompareParameters bondCompParams;
        bondCompParams.MatchStereo = true;
        bondCompParams.RingMatchesRingOnly = false;
        bondCompParams.CompleteRingsOnly = true;
        bondCompParams.MatchFusedRings = true;
        bondCompParams.MatchFusedRingsStrict = false;
        mcsParams.BondCompareParameters = bondCompParams;
        mcsParams.Timeout = 1;
        mcsParams.setMCSAtomTyperFromEnum(RDKit::AtomCompareElements);
        mcsParams.setMCSBondTyperFromEnum(RDKit::BondCompareOrderExact);
        return mcsParams;
    }
}  // namespace coaler::core