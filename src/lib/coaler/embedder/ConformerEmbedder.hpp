//
// Created by chris on 12/4/23.
//

#pragma once
#include <GraphMol/ROMol.h>

namespace coaler::embedder{

    using CoreAtomMapping = std::map<int, RDGeom::Point3D>;

    class ConformerEmbedder {
      public:

        ConformerEmbedder(const RDKit::ROMol& core,
                          unsigned numCoreConfigs);

        ConformerEmbedder(const RDKit::ROMol& core,
                          const CoreAtomMapping& coreMap);

        void embedWithFixedCore(
            const RDKit::ROMol& mol,
            const RDKit::ROMol& core
        );

      private:
        RDKit::ROMol m_core;
    };
}

