/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#pragma once
#include <GraphMol/ROMol.h>

namespace coaler::embedder {

    using CoreAtomMapping = std::map<int, RDGeom::Point3D>;

    /**
     * The ConformerEmbedder class provides functionality for the generation of conformers for
     * a given molecule with contrained core coordinates.
     */
    class ConformerEmbedder {
      public:
        ConformerEmbedder(RDKit::ROMol core, const CoreAtomMapping& coreMap);

        /**
         * Embed a fixed number of conformes at the first core match.
         * @return True upon success.
         */
        bool embedWithFixedCore(RDKit::ROMol& mol, unsigned numConfs);

        /**
         * Embed an even amount of Conformers at every core match.
         * @param mol The molecule to embed.
         *
         * @note If the core has a too high symmetry, it is possible, that no embedding can be
         * performed within the given min/max constraints.
         *
         * @return True upon success.
         */
        bool embedEvenlyAcrossAllMatches(RDKit::ROMol& mol, unsigned minNofConfs, unsigned maxNofConfs);

        static std::vector<unsigned> distributeApproxEvenly(unsigned nofMatches, unsigned maxConformers);

      private:
        RDKit::ROMol m_core;
    };
}  // namespace coaler::embedder
