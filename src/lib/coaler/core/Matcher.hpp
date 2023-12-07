//
// Created by malte on 12/4/23.
//

#pragma once

#include "Forward.hpp"

namespace coaler::core {

    class Matcher {
      public:
        /**
         * calculates the MCS of the molecules
         * @return MCS as ROMol
         */
        static std::optional<RDKit::ROMOL_SPTR> calculateCoreMcs(RDKit::MOL_SPTR_VECT mols);

        /**
         * calculates the Murcko scaffold of the molecules
         * @return Murcko Scaffold as ROMol
         */
        static std::optional<RDKit::ROMOL_SPTR> calculateCoreMurcko(RDKit::MOL_SPTR_VECT mols);
    };
}  // namespace coaler::core