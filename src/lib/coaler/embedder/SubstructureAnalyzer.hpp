//
// Created by chris on 12/7/23.
//

#pragma once

#include <GraphMol/ROMol.h>

namespace coaler::embedder {

    class SubstructureAnalyzer {
      public:
        static unsigned getNumberOfRingRotations(const RDKit::ROMol& molecule);

        static unsigned getNumberOfUniqueSubstructureMatches(const RDKit::ROMol& query, const RDKit::ROMol& molecule);
    };
}  // namespace coaler::embedder
