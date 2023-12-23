//
// Created by malte on 12/4/23.
//

#include "Matcher.hpp"

#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <spdlog/spdlog.h>

#include "GraphMol/MolOps.h"
#include "GraphMol/ChemTransforms/ChemTransforms.h"
#include "GraphMol/DistGeomHelpers/Embedder.h"
#include "GraphMol/FMCS/FMCS.h"
#include "unordered_set"

namespace coaler::core {
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


        RDKit::SubstructMatchParameters substructMatchParams;
        // setting this to true breaks SIENA 2w0v
        substructMatchParams.useChirality = false;
        substructMatchParams.useEnhancedStereo = true;
        substructMatchParams.aromaticMatchesConjugated = true;
        substructMatchParams.numThreads = numOfThreads;

        std::vector<RDKit::MatchVectType> structMatches
            = RDKit::SubstructMatch(first, *mcs.QueryMol, substructMatchParams);


        assert(!structMatches.empty());

        std::unordered_set<uint> matchedAtoms;
        for (auto const &match: structMatches) {
            for (auto const &[queryId, molId]: match) {
                matchedAtoms.insert(molId);
            }
        }

        std::vector<uint> atomsForRemoval;
        for (auto i = first.getNumAtoms()-1; i == 0; i--) {
            if (!matchedAtoms.count(i)) {
                atomsForRemoval.push_back(i);
            }
        }

        for (auto const &atomId: atomsForRemoval) {
            first.removeAtom(atomId);
        }

        auto params = RDKit::DGeomHelpers::srETKDGv3;
        params.numThreads = numOfThreads;
        RDKit::DGeomHelpers::EmbedMolecule(first, params);

        std::vector<std::pair<int, double>> result;
        RDKit::MMFF::MMFFOptimizeMoleculeConfs(first,result,numOfThreads);

        return CoreResult{mcs.QueryMol, boost::make_shared<RDKit::ROMol>(first)};
    }

    std::optional<CoreResult> Matcher::calculateCoreMurcko(RDKit::MOL_SPTR_VECT mols, int numOfThreads) {
        auto mcs = Matcher::calculateCoreMcs(mols, numOfThreads);
        if (!mcs.has_value()) {
            return std::nullopt;
        }

        RDKit::ROMol ret = *RDKit::MurckoDecompose(*mols.at(0));
        return CoreResult{};
    }

}  // namespace coaler::core