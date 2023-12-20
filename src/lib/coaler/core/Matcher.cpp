//
// Created by malte on 12/4/23.
//

#include "Matcher.hpp"

namespace coaler::core {
//    std::optional<CoreResult> Matcher::calculateCoreMurcko(RDKit::MOL_SPTR_VECT mols) {
//        auto mcs = Matcher::calculateCoreMcs(mols);
//        if (!mcs.has_value()) {
//            return std::nullopt;
//        }
//
//        RDKit::ROMol ret = *RDKit::MurckoDecompose(*mols.at(0));
//        return std::make_pair(nullptr, mcs.value().second);
//    }

    std::optional<CoreResult> Matcher::getResult() {
        return this->core_;
    }

    std::optional<CoreMatch> Matcher::matchMolecule(RDKit::ROMOL_SPTR mol) {
        const RDKit::SubstructMatchParameters matchParams;
        auto matches = RDKit::SubstructMatch(*mol, *this->core_.first, matchParams);
        if (matches.empty()) {
            return std::nullopt;
        }

        return matches;
    }
}  // namespace coaler::core