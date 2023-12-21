//
// Created by chris on 12/21/23.
//

#pragma once
#include <GraphMol/DistGeomHelpers/Embedder.h>

#include <map>

namespace coaler::embedder {

    struct ConformerEmbeddingParams {
        unsigned minConfsPerMatch;
        unsigned maxConfsPerMatch;
        unsigned maxTotalConfsPerMol;
    };

    using CoreAtomMapping = std::map<int, RDGeom::Point3D>;

}  // namespace coaler::embedder