#pragma once

#include <GraphMol/ROMol.h>

namespace coaler::embedder {

    class [[maybe_unused]] CoreSymmetryCalculator {
      public:
        [[maybe_unused]] static unsigned getNofSymmetryAxes(const RDKit::ROMol& mol);
    };

}  // namespace coaler::embedder
