//
// Created by malte on 12/4/23.
//

#include "Matcher.hpp"

#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>

#include "GraphMol/ChemTransforms/ChemTransforms.h"
#include "GraphMol/DistGeomHelpers/Embedder.h"
#include "GraphMol/FMCS/FMCS.h"

namespace coaler::core {
    struct CompareUnsignedInts {
        bool operator()(unsigned int a, unsigned int b) const { return a > b; }
    };

    std::optional<RDKit::ROMOL_SPTR> Matcher::calculateCoreMcs(RDKit::MOL_SPTR_VECT mols) {
        RDKit::MCSParameters mcsParams;
        RDKit::MCSAtomCompareParameters atomCompParams;
        atomCompParams.MatchChiralTag = true;
        atomCompParams.MatchFormalCharge = true;
        atomCompParams.MatchIsotope = true;
        atomCompParams.MatchValences = true;
        atomCompParams.RingMatchesRingOnly = true;
        atomCompParams.CompleteRingsOnly = true;
        mcsParams.AtomCompareParameters = atomCompParams;

        RDKit::MCSBondCompareParameters bondCompParams;
        bondCompParams.MatchStereo = true;
        bondCompParams.RingMatchesRingOnly = true;
        bondCompParams.CompleteRingsOnly = true;
        bondCompParams.MatchFusedRingsStrict = true;
        mcsParams.BondCompareParameters = bondCompParams;

        RDKit::MCSResult const mcs = RDKit::findMCS(mols, &mcsParams);
        if (mcs.QueryMol == nullptr) {
            return std::nullopt;
        }

        RDKit::RWMol first = *mols.at(0);

        auto structMatches = RDKit::SubstructMatch(first, *mcs.QueryMol);
        assert(structMatches.size() > 0);

        std::set<unsigned> atomIds;
        for (auto i = 0; i < first.getNumAtoms(); i++) {
            atomIds.emplace(i);
        }

        std::set<unsigned> match;
        // only search first match
        for (auto matchPair : structMatches.at(0)) {
            match.emplace(matchPair.second);
        }

        std::set<unsigned, CompareUnsignedInts> difference;
        std::set_difference(atomIds.begin(), atomIds.end(), match.begin(), match.end(),
                            std::inserter(difference, difference.begin()));

        for (auto id : difference) {
            first.removeAtom(id);
        }

        RDKit::DGeomHelpers::EmbedParameters params;
        RDKit::DGeomHelpers::EmbedMolecule(first, params);

        RDKit::MolOps::sanitizeMol(first);

        return boost::make_shared<RDKit::ROMol>(first);
    }

    std::optional<RDKit::ROMOL_SPTR> Matcher::calculateCoreMurcko(RDKit::MOL_SPTR_VECT mols) {
        auto mcs = Matcher::calculateCoreMcs(mols);
        if (!mcs.has_value()) {
            return std::nullopt;
        }

        RDKit::ROMol ret = *RDKit::MurckoDecompose(*mcs.value());
        return boost::make_shared<RDKit::ROMol>(ret);
    }

}  // namespace coaler::core