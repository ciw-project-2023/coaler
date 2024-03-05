#pragma once

#include <GraphMol/ROMol.h>

namespace coaler::embedder {

    /**
     * The SubstructureAnalyzer class provides functionality for the analysis of substructures in molecules.
     */
    class SubstructureAnalyzer {
      public:
        /**
         * Calculate the number of identity core rotations.
         * c1ccccc1 --> 6
         * c1cnccn1 --> 2
         *
         * @note This function only works on single rings. If not all atoms have degree 2, the
         * default (1) is returned.
         *
         * @param Molecule the ring molecule to analyze.
         * @return The number of identity core rotations.
         */
        static unsigned getNumberOfRingRotations(const RDKit::ROMol& molecule);

        static unsigned getNumberOfUniqueSubstructureMatches(const RDKit::ROMol& query, const RDKit::ROMol& molecule);
    };
}  // namespace coaler::embedder
