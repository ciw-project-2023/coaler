//
// Created by chris on 12/16/23.
//

#include "BucketDistributor.hpp"
std::vector<unsigned int> coaler::embedder::BucketDistributor::distributeApproxEvenly(unsigned int nofMatches,
                                                                                      unsigned int maxConformers) {
    std::vector<unsigned> confsForMatch(nofMatches);
    unsigned decrementPosition = maxConformers % nofMatches;
    unsigned baseNofConfs = maxConformers / nofMatches;

    std::fill(confsForMatch.begin(), confsForMatch.begin() + decrementPosition, baseNofConfs + 1);

    std::fill(confsForMatch.begin() + decrementPosition, confsForMatch.end(), baseNofConfs);

    return confsForMatch;
}
