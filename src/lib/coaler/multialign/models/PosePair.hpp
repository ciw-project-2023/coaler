#pragma once
#include "Alias.hpp"
#include "UniquePoseID.hpp"
#include "boost/functional/hash.hpp"
/*!
 * @file PosePair.hpp
 * @brief This file contains the PosePair class which is used to represent a pair of poses.
 */
namespace coaler::multialign {

    /**
     * @brief PosePair class to represent a pair of poses.
     *
     * The PosePair class is used to represent a pair of poses.
     * It is used to store the poses in a single object.
     */
    class PosePair {
      public:
        /**
         * @brief Constructor for the PosePair class.
         *
         * @param first The first pose.
         * @param second The second pose.
         */
        explicit PosePair(UniquePoseID first, UniquePoseID second);

        [[nodiscard]] UniquePoseID getFirst() const noexcept;
        [[nodiscard]] UniquePoseID getSecond() const noexcept;

        bool operator==(const PosePair& other) const;

      private:
        UniquePoseID m_firstPose;
        UniquePoseID m_secondPose;
    };

    /**
     * @brief Hash function for PosePair
     */
    struct PosePairHash {
        std::size_t operator()(const PosePair& pair) const {
            std::size_t seed = 0;
            boost::hash_combine(seed, UniquePoseIdentifierHash()(pair.getFirst()));
            boost::hash_combine(seed, UniquePoseIdentifierHash()(pair.getSecond()));
            return seed;
        }
    };

}  // namespace coaler::multialign
