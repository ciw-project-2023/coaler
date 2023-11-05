//
// Created by chris on 11/5/23.
//

#pragma once
#include "Forward.hpp"

namespace MultiAlign
{

    class PosePair{
    public:
        PosePair(PoseID first,
                 PoseID second);

        [[nodiscard]] const PoseID getFirst() const noexcept;
        [[nodiscard]] const PoseID getSecond() const noexcept;

        bool operator==(const PosePair& other) const;

    private:
        PoseID m_firstPose;
        PoseID m_secondPose;
    };

    struct PosePairHash
    {
        std::size_t operator()(const PosePair& pair) const
        {
            std::size_t seed = 0;

            boost::hash_combine(seed,boost::hash_value(pair.getFirst()));
            boost::hash_combine(seed,boost::hash_value(pair.getSecond()));

            return seed;
        }};

}


