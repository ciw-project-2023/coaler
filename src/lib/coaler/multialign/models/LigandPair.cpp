/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#include "LigandPair.hpp"

#include <cassert>

namespace coaler::multialign {
    LigandPair::LigandPair(const LigandID first, const LigandID second) {
        // swap values to rule out duplicates
        assert(first != second);
        if (first > second) {
            this->first = second;
            this->second = first;
        } else {
            this->first = first;
            this->second = second;
        }
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    bool LigandPair::operator==(const LigandPair &other) const {
        return this->first == other.first && this->second == other.second;
    }
}  // namespace coaler::multialign
