/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#include "PosePair.hpp"

#include <boost/functional/hash.hpp>
#include <cassert>

namespace coaler::multialign {
    PosePair::PosePair(const UniquePoseID first, const UniquePoseID second) {
        // swap values to rule out duplicates
        assert(first != second);
        if (first > second) {
            m_firstPose = second;
            m_secondPose = first;
        } else {
            m_firstPose = first;
            m_secondPose = second;
        }
    }

    bool PosePair::operator==(const PosePair &other) const {
        return this->m_firstPose == other.m_firstPose && this->m_secondPose == other.m_secondPose;
    }

    const UniquePoseID PosePair::getFirst() const noexcept { return m_firstPose; }

    const UniquePoseID PosePair::getSecond() const noexcept { return m_secondPose; }
}  // namespace coaler::multialign
