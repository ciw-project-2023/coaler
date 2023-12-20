//
// Created by niklas on 12/18/23.
//

#include "MatcherMCS.h"

namespace coaler::core {
    std::optional<RDKit::MCSResult> MatcherMCS::findMCS(RDKit::MOL_SPTR_VECT mols) {
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

        RDKit::MCSResult mcs = RDKit::findMCS(mols, &mcsParams);
        if (mcs.QueryMol == nullptr) {
            return std::nullopt;
        }

        return mcs;
    }

    MatcherMCS::MatcherMCS(RDKit::MOL_SPTR_VECT mols) : Matcher() {
        auto mcsResult = this->findMCS(mols);
        if (!mcsResult.has_value()) {
            throw std::runtime_error("Could not find MCS");
        }

        auto mcs = mcsResult.value();

        spdlog::info("MCS: {}", mcs.SmartsString);
        RDKit::RWMol first = *mols.at(0);

        auto params = RDKit::DGeomHelpers::srETKDGv3;
        params.randomSeed = 42;
        params.useRandomCoords = true;
        RDKit::DGeomHelpers::EmbedMolecule(first, params);

        auto structMatches = RDKit::SubstructMatch(first, *mcs.QueryMol);
        assert(!structMatches.empty());

        AtomMap moleculeCoreCoords;
        RDKit::Conformer conformer = first.getConformer(0);
        for (const auto& [queryId, molId] : structMatches.at(0)) {
            const RDGeom::Point3D atomCoords = conformer.getAtomPos(molId);
            moleculeCoreCoords.emplace(queryId, atomCoords);
        }

        this->core_ = std::make_pair(mcs.QueryMol, moleculeCoreCoords);
    }
}
