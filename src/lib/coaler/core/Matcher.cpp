//
// Created by malte on 12/4/23.
//

#include "Matcher.hpp"

#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <spdlog/spdlog.h>

#include "GraphMol/ChemTransforms/ChemTransforms.h"
#include "GraphMol/DistGeomHelpers/Embedder.h"
#include "GraphMol/FMCS/FMCS.h"

namespace coaler::core {
    std::optional<CoreResult> Matcher::calculateCoreMcs(RDKit::MOL_SPTR_VECT mols, int numOfThreads) {
        // Generates all parameters needed for RDKit::findMCS()
        const RDKit::MCSParameters mcsParams = getMCSParams();

        RDKit::MCSResult const mcs = RDKit::findMCS(mols, &mcsParams);
        if (mcs.QueryMol == nullptr) {
            return std::nullopt;
        }

        spdlog::info("MCS: {}", mcs.SmartsString);

        RDKit::RWMol first = *mols.at(0);

        auto params = RDKit::DGeomHelpers::srETKDGv3;
        RDKit::DGeomHelpers::EmbedMolecule(first, params);

        std::vector<std::pair<int, double>> result;
        RDKit::MMFF::MMFFOptimizeMoleculeConfs(first, result);

        const RDKit::SubstructMatchParameters substructMatchParams = getSubstructMatchParams(numOfThreads);

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

        RDKit::ROMol ret = *RDKit::MurckoDecompose(*mols.at(0));
        auto ret_ptr = boost::make_shared<RDKit::ROMol>(ret);
        return std::make_pair(ret_ptr, mcs.value().second);
    }

    RDKit::MCSParameters Matcher::getMCSParams() {
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
        return mcsParams;
    }

    RDKit::SubstructMatchParameters Matcher::getSubstructMatchParams(int numOfThreads) {
        RDKit::SubstructMatchParameters substructMatchParams;
        substructMatchParams.uniquify = true;
        substructMatchParams.useChirality = true;
        substructMatchParams.useEnhancedStereo = true;
        substructMatchParams.aromaticMatchesConjugated = true;
        substructMatchParams.numThreads = numOfThreads;
        return substructMatchParams;
    }

}  // namespace coaler::core