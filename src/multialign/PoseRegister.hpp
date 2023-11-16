//
// Created by chris on 11/4/23.
//

#pragma once

#include "Forward.hpp"
#include "BasicClasses/PosePair.hpp"
#include <unordered_map>

namespace MultiAlign
{
    using RegisterMap = std::unordered_map<PosePair, double, PosePairHash>;

    class PoseRegister {

    public:

        PoseRegister(LigandID firstLigand,
                              LigandID secondLigand,
                              unsigned maxSize);

        [[nodiscard]] unsigned getSize() const noexcept;

        void setMaxSize(unsigned max) noexcept;

        [[nodiscard]] LigandID getFirstLigandID() const noexcept;

        [[nodiscard]] LigandID getSecondLigandID() const noexcept;

        bool addPoses(PosePair pair, double score);

        PosePair getMinimumPosePairInRegister();

    private:

        void updateLowestScoringPosePair();

        LigandID m_first;
        LigandID m_second;
        unsigned m_maxSize;
        RegisterMap m_register;
        PosePair m_lowestScoringPosePair;
    };

}

