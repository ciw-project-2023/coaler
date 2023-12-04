/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */
#pragma once
#include "Forward.hpp"
#include "UniquePoseID.hpp"

namespace coaler::multialign {

    /**
     * A pair of conformers.
     */
    class PosePair {
      public:
        explicit PosePair(UniquePoseID first, UniquePoseID second);

        [[nodiscard]] const UniquePoseID getFirst() const noexcept;
        [[nodiscard]] const UniquePoseID getSecond() const noexcept;

        bool operator==(const PosePair& other) const;

      private:
        UniquePoseID m_firstPose;
        UniquePoseID m_secondPose;
    };

    struct PosePairHash {
        std::size_t operator()(const PosePair& pair) const {
            std::size_t seed = 0;
            boost::hash_combine(seed, UniquePoseIdentifierHash()(pair.getFirst()));
            boost::hash_combine(seed, UniquePoseIdentifierHash()(pair.getSecond()));
            return seed;
        }
    };

}  // namespace coaler::multialign
