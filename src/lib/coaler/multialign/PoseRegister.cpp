/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#include "PoseRegister.hpp"

#include "models/Forward.hpp"

namespace {
    struct PosePairScoreGreater {
        bool operator()(const PosePairAndScore lhs, const PosePairAndScore rhs) { return lhs.second > rhs.second; }
    };

}  // namespace

namespace coaler::multialign {

    PoseRegister::PoseRegister(LigandID firstLigand, LigandID secondLigand, unsigned maxSize)
        : m_maxSize(maxSize), m_first(firstLigand), m_second(secondLigand),
          m_lowest(std::make_pair(PosePair({0,1},{0,0}), std::numeric_limits<double>::max())),
          m_highest(std::make_pair(PosePair({0,1},{0,0}), std::numeric_limits<double>::min()))
    {
        assert(m_first != m_second);
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    unsigned PoseRegister::getSize() const noexcept { return m_register.size(); }

    LigandID PoseRegister::getFirstLigandID() const noexcept { return m_first; }

    LigandID PoseRegister::getSecondLigandID() const noexcept { return m_second; }

    /*----------------------------------------------------------------------------------------------------------------*/

    void PoseRegister::addPoses(const PosePair pair, const double score) {
        if (m_register.size() < m_maxSize) {
            m_register.emplace_back(pair, score);
            PosePairAndScore insertedElement = m_register.back();
            updateHighest(insertedElement);
            updateLowest();
        } else if (score > m_lowest.second) {
            m_register.emplace_back(pair, score);
            updateHighest(m_register.back());
            removeAndUpdateLowest();
        }
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    PosePair PoseRegister::getHighestScoringPair() { return m_highest.first; }

    /*----------------------------------------------------------------------------------------------------------------*/

    bool PoseRegister::containsPose(const UniquePoseID &pose) {
        return std::any_of(m_register.begin(), m_register.end(), [pose](const auto entry) {
            return entry.first.getFirst() == pose || entry.first.getSecond() == pose;
        });
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    PosePair PoseRegister::getHighestScoringPosePairForPose(const UniquePoseID &pose) {
        assert(this->containsPose(pose));
        PosePairAndScore currentBest = std::make_pair(m_register.at(0).first, -1);  // dummy value, will be overwritten
        for (const PosePairAndScore &entry : m_register) {
            if (entry.first.getFirst() == pose || entry.first.getSecond() == pose) {
                if (entry.second > currentBest.second) {
                    currentBest = entry;
                }
            }
        }
        return currentBest.first;
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    void PoseRegister::updateLowest() {
        auto lowest = std::min_element(m_register.begin(), m_register.end(), PosePairScoreGreater());
        m_lowest = *lowest;
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    void PoseRegister::updateHighest(const PosePairAndScore& insertedPair) {
        if (insertedPair.second > m_highest.second) {
            m_highest = insertedPair;
        }
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    void PoseRegister::removeAndUpdateLowest() {
        assert(m_register.size() == m_maxSize + 1);
        auto itemToRemove = m_register.begin();
        auto newLowest = m_register.begin();
        for(auto iter = m_register.begin(); iter != m_register.end(); iter++) {
            if(iter->second < itemToRemove->second){
                newLowest = itemToRemove;
                itemToRemove = iter;
            }
            else if(iter->second < newLowest->second) {
                newLowest = iter;
            }
        }
        m_register.erase(itemToRemove);
        m_lowest = *newLowest;
        assert(m_register.size() == m_maxSize);

    }

}  // namespace coaler::multialign
