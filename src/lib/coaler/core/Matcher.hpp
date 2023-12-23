//
// Created by malte on 12/4/23.
//

#pragma once

#include "Forward.hpp"

namespace coaler::core {
    struct CoreResult {
        RDKit::ROMOL_SPTR core;
        RDKit::ROMOL_SPTR ref;
    };

    class Matcher {
      public:
        /**
         * calculates the MCS of the molecules
         * @return MCS as ROMol
         */
        static std::optional<CoreResult> calculateCoreMcs(RDKit::MOL_SPTR_VECT mols, int numOfThreads);

        /**
         * calculates the Murcko scaffold of the molecules
         * @return Murcko Scaffold as ROMol
         */
        static std::optional<CoreResult> calculateCoreMurcko(RDKit::MOL_SPTR_VECT mols, int numOfThreads);
    };
}  // namespace coaler::core