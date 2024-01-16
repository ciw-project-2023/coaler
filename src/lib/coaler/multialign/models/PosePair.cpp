#include "PosePair.hpp"

#include <boost/functional/hash.hpp>
#include <cassert>

namespace coaler::multialign {
    // NOLINTBEGIN(cppcoreguidelines-pro-type-member-init)
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
    // NOLINTEND(cppcoreguidelines-pro-type-member-init)

    bool PosePair::operator==(const PosePair &other) const {
        return this->m_firstPose == other.m_firstPose && this->m_secondPose == other.m_secondPose;
    }

    UniquePoseID PosePair::getFirst() const noexcept { return m_firstPose; }

    UniquePoseID PosePair::getSecond() const noexcept { return m_secondPose; }
}  // namespace coaler::multialign
