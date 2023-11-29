//
// Created by chris on 11/5/23.
//

#pragma once
#include "Forward.hpp"
#include "UniquePoseIdentifier.hpp"

namespace coaler::multialign
{

    /**
     * A pair of conformers.
     */
    class PosePair{
    public:
        explicit PosePair(UniquePoseIdentifier first,
                 UniquePoseIdentifier second);

        [[nodiscard]] const UniquePoseIdentifier getFirst() const noexcept;
        [[nodiscard]] const UniquePoseIdentifier getSecond() const noexcept;

        bool operator==(const PosePair& other) const;

    private:
        UniquePoseIdentifier m_firstPose;
        UniquePoseIdentifier m_secondPose;
    };

    struct PosePairHash
    {
        std::size_t operator()(const PosePair& pair) const
        {
            std::size_t seed = 0;
            boost::hash_combine(seed, UniquePoseIdentifierHash()(pair.getFirst()));
            boost::hash_combine(seed, UniquePoseIdentifierHash()(pair.getSecond()));
            return seed;
        }};

}


