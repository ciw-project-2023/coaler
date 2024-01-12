#pragma once

#include <GraphMol/ROMol.h>

namespace coaler::embedder {
    class CoreSymmetryCalculator {
      public:
        static unsigned getNofSymmetryAxes(const RDKit::ROMol& mol);
    };
}  // namespace coaler::embedder
