//
// Created by chris on 12/21/23.
//

#pragma once

namespace coaler::embedder{
    struct ConformerEmbeddingParams{
        unsigned maxTotalConfsPerMol;
        unsigned minConfsPerMatch;
        unsigned maxConfsPerMatch;
    };
}