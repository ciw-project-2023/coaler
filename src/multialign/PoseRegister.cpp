//
// Created by chris on 11/4/23.
//

#include "PoseRegister.hpp"

#include "BasicClasses/Forward.hpp"

namespace {
    struct PosePairScoreGreater {
        bool operator()(const PosePairAndScore lhs, const PosePairAndScore rhs) { return lhs.second > rhs.second; }
    };

}  // namespace

namespace coaler::multialign {

    PoseRegister::PoseRegister(LigandID firstLigand, LigandID secondLigand, unsigned maxSize)
        : m_maxSize(maxSize), m_first(firstLigand), m_second(secondLigand) {
        assert(m_first != m_second);
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    unsigned PoseRegister::getSize() const noexcept { return m_register.size(); }

    LigandID PoseRegister::getFirstLigandID() const noexcept { return m_first; }

    LigandID PoseRegister::getSecondLigandID() const noexcept { return m_second; }

    /*----------------------------------------------------------------------------------------------------------------*/

    bool PoseRegister::addPoses(const PosePair pair, const double score) {
        if (m_register.size() < m_maxSize) {
            m_register.emplace_back(pair, score);
            std::sort(m_register.begin(), m_register.end(), PosePairScoreGreater());
            return true;
        } else if (score > m_register.back().second) {
            m_register.pop_back();
            m_register.emplace_back(pair, score);
            std::sort(m_register.begin(), m_register.end(), PosePairScoreGreater());
        }
        return false;
    }

    PosePair PoseRegister::getHighestScoringPair() { return m_register.at(0).first; }

    /*----------------------------------------------------------------------------------------------------------------*/

    bool PoseRegister::containsPose(const UniquePoseIdentifier &pose) {
        return std::any_of(m_register.begin(), m_register.end(), [pose](const auto entry) {
            return entry.first.getFirst() == pose || entry.first.getSecond() == pose;
        });
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    PosePair PoseRegister::getHighestScoringPosePairForPose(const UniquePoseIdentifier &pose) {
        assert(this->containsPose(pose));
        PosePairAndScore currentBest = std::make_pair(m_register.at(0).first, -1); //dummy value, will be overwritten
        for(const PosePairAndScore& entry : m_register)
        {
            if(entry.first.getFirst() == pose || entry.first.getSecond() == pose)
            {
                if(entry.second > currentBest.second)
                {
                    currentBest = entry;
                }
            }
        }
        return currentBest.first;
    }

}  // namespace coaler::multialign