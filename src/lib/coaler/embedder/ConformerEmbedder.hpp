/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#pragma once
#include <GraphMol/ROMol.h>

namespace coaler::embedder {

    using CoreAtomMapping = std::map<int, RDGeom::Point3D>;

    class ConformerEmbedder {
      public:

        explicit ConformerEmbedder(RDKit::ROMOL_SPTR& core, int threads = 1);

        void embedConformersWithFixedCore(RDKit::ROMOL_SPTR mol, unsigned numConfs);

        bool embedEvenlyAcrossAllMatches(RDKit::ROMol& mol, unsigned minNofConfs, unsigned maxNofConfs);

        static std::vector<unsigned> distributeApproxEvenly(unsigned nofMatches, unsigned maxConformers);

      private:

        RDKit::ROMOL_SPTR m_core;
        int m_threads;

    };
}  // namespace coaler::embedder
