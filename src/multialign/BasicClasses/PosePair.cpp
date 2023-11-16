//
// Created by chris on 11/5/23.
//

#include "PosePair.hpp"

namespace MultiAlign
{
    PosePair::PosePair(const PoseID first, const PoseID second) {
        //swap values to rule out duplicates
        assert(first != second);
        if(first > second){
            m_firstPose = second;
            m_secondPose = first;
        } else {
            m_firstPose = first;
            m_secondPose = second;
        }
    }

    bool PosePair::operator==(const PosePair &other) const {
        return this->m_firstPose == other.m_firstPose
        && this->m_secondPose == other.m_secondPose;
    }

    const PoseID PosePair::getFirst() const noexcept {
        return m_firstPose;
    }

    const PoseID PosePair::getSecond() const noexcept {
        return m_secondPose;
    }
}