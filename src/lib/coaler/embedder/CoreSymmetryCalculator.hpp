#pragma once

#include <GraphMol/ROMol.h>

#include "Forward.hpp"

namespace coaler::embedder {

    class CoreSymmetryCalculator {
      public:
        static unsigned getNofSymmetryAxes(const RDKit::ROMol& mol);

        static CoreAtomMapping getShiftedMapping(const coaler::embedder::CoreAtomMapping& map, unsigned shift);

        static std::vector<unsigned> distributeApproxEvenly(unsigned int nofMatches, unsigned int maxConformers,unsigned int maxPerMatch);
    };
} // namespace coaler::embedder
