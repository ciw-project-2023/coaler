/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */
#pragma once
#include <cstddef>

#include "Alias.hpp"
#include "boost/functional/hash.hpp"

namespace coaler::multialign {
    /**
     * A pair of ligands.
     */
    class LigandPair : public std::pair<LigandID, LigandID> {
      public:
        using std::pair<LigandID, LigandID>::pair;

        explicit LigandPair(LigandID first, LigandID second);

        bool operator==(const LigandPair& other) const;
    };

    struct LigandPairHash {
        std::size_t operator()(const LigandPair& pair) const {
            std::size_t seed = 0;

            boost::hash_combine(seed, boost::hash_value(std::get<0>(pair)));
            boost::hash_combine(seed, boost::hash_value(std::get<1>(pair)));

            return seed;
        }
    };
}  // namespace coaler::multialign
