/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#pragma once
#include <GraphMol/ROMol.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>

namespace coaler::embedder {

    using CoreAtomMapping = std::map<int, RDGeom::Point3D>;

    /**
     * The ConformerEmbedder class provides functionality for the generation of conformers for
     * a given molecule with contrained core coordinates.
     */
    class ConformerEmbedder {
      public:
        explicit ConformerEmbedder(RDKit::ROMOL_SPTR& query, CoreAtomMapping& coords, int threads = 1, bool divideConformersByMatches = false);

        /**
         * Embed an even amount of Conformers at every core match.
         * @param mol The molecule to embed.
         *
         * @note If the core has a too high symmetry, it is possible, that no embedding can be
         * performed within the given min/max constraints.
         *
         * @return True upon success.
         */
        void embedEvenlyAcrossAllMatches(const RDKit::ROMOL_SPTR &mol, unsigned numConfs);

      private:
        RDKit::ROMOL_SPTR m_core;
        CoreAtomMapping m_coords;
        int m_threads;
        bool m_divideConformersByMatches;


        RDKit::DGeomHelpers::EmbedParameters getEmbeddingParameters(const CoreAtomMapping& coords);
    };
}  // namespace coaler::embedder
