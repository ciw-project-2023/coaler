//
// Created by chris on 12/4/23.
//

#include "ConformerEmbedder.hpp"

#include <GraphMol/DistGeomHelpers/Embedder.h>

namespace coaler::embedder
{
    ConformerEmbedder::ConformerEmbedder(const RDKit::ROMol& core, const coaler::embedder::CoreAtomMapping& coreMap)
    : m_core(core){
        RDKit::DGeomHelpers::EmbedParameters params;
        params.randomSeed = 42;
        params.coordMap = &coreMap;
    }
}

