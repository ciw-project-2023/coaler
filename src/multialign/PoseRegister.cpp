//
// Created by chris on 11/4/23.
//

#include "PoseRegister.hpp"
#include "BasicClasses/Forward.hpp"

namespace MultiAlign {


    PoseRegister::PoseRegister(LigandID firstLigand,
                               LigandID secondLigand,
                               unsigned maxSize)
    : m_maxSize(maxSize)
    , m_first(firstLigand)
    , m_second(secondLigand)
    , m_lowestScoringPosePair(
            std::numeric_limits<unsigned>::max(),
            std::numeric_limits<unsigned>::max() - 1)
    {
        assert(m_first != m_second);
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    unsigned PoseRegister::getSize() const noexcept {
        return m_register.size();
    }

    LigandID PoseRegister::getFirstLigandID() const noexcept {
        return m_first;
    }

    LigandID PoseRegister::getSecondLigandID() const noexcept {
        return m_second;
    }

    void PoseRegister::setMaxSize(unsigned int max) noexcept {
        m_maxSize = max;
    }

    bool PoseRegister::addPoses(const PosePair pair, const double score) {
        if(m_register.size() < m_maxSize)
        {
            m_register.insert({pair, score});
            updateLowestScoringPosePair(); //NOTE can also check directly
            return true;
        }
        else if(score > m_register.at(m_lowestScoringPosePair))
        {
            m_register.erase(m_lowestScoringPosePair);
            m_register.insert({pair, score});
            updateLowestScoringPosePair();
            return true;
        }
        return false;
    }

    void PoseRegister::updateLowestScoringPosePair() {
        // if only one pose pair in register its automatically the lowest
        if(m_register.size() == 1)
        {
            m_lowestScoringPosePair = m_register.begin()->first;
            return;
        }

        //else find min pair
        auto min = m_register.begin();
        for(auto iter = m_register.begin(); iter != m_register.end(); iter++)
        {
            if(iter->second < min->second)
            {
                min = iter;
            }
        }
        m_lowestScoringPosePair = min->first;
    }

    PosePair PoseRegister::getMinimumPosePairInRegister() {
        return m_lowestScoringPosePair;
    }

}