//
// Created by malte on 12/4/23.
//

#pragma once

#include "Forward.hpp"

namespace coaler::core {
    using AtomMap = std::map<int, RDGeom::Point3D>;
    using CoreResult = std::pair<RDKit::ROMOL_SPTR, AtomMap>;

    class Matcher {
      public:
        /**
         * calculates the MCS of the molecules
         * @return MCS as ROMol
         */
        static std::optional<CoreResult> calculateCoreMcs(RDKit::MOL_SPTR_VECT mols, int numOfThreads = 1);

        /**
         * calculates the Murcko scaffold of the molecules
         * @return Murcko Scaffold as ROMol
         */
        static std::optional<CoreResult> calculateCoreMurcko(RDKit::MOL_SPTR_VECT mols, int numOfThreads = 1);
    };
}  // namespace coaler::core