//
// Created by chris on 12/16/23.
//

#pragma once

#include <vector>
namespace coaler::embedder {

    class BucketDistributor {
      public:
        static std::vector<unsigned> distributeApproxEvenly(unsigned nofMatches, unsigned maxConformers);
    };
}  // namespace coaler::embedder
