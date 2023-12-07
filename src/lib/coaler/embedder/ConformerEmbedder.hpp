/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#pragma once
#include <GraphMol/ROMol.h>

namespace coaler::embedder {

    using CoreAtomMapping = std::map<int, RDGeom::Point3D>;

    class ConformerEmbedder {
      public:

        ConformerEmbedder(RDKit::ROMol core, const CoreAtomMapping& coreMap);

        bool embedWithFixedCore(RDKit::ROMol& mol, unsigned numConfs);

        bool embedEvenlyAcrossAllMatches(RDKit::ROMol& mol, unsigned minNofConfs, unsigned maxNofConfs);

        static std::vector<unsigned> distributeApproxEvenly(unsigned nofMatches, unsigned maxConformers);

      private:

        RDKit::ROMol m_core;
    };
}  // namespace coaler::embedder
