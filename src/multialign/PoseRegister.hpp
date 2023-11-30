//
// Created by chris on 11/4/23.
//

#pragma once

#include <unordered_map>

#include "BasicClasses/PosePair.hpp"
#include "Forward.hpp"

using PosePairAndScore = std::pair<coaler::multialign::PosePair, double>;

namespace coaler::multialign {

    /**
     * For a pair of ligands this contains the best aligning pairwise poses.
     */
    class PoseRegister {
      public:
        PoseRegister(LigandID firstLigand, LigandID secondLigand, unsigned maxSize);

        /**
         *
         * @return The number of aligned pose pairs in the register.
         */
        [[nodiscard]] unsigned getSize() const noexcept;

        [[nodiscard]] LigandID getFirstLigandID() const noexcept;

        [[nodiscard]] LigandID getSecondLigandID() const noexcept;

        /**
         * Adds a pose pair to the register. If the register is full
         * and the pose is inferior to the worst scoring pose in the register,
         * it is discarded.
         * @param pair The pose pair to add.
         * @param score The score of the pairs alignment.
         * @return True if the pose was added.
         */
        bool addPoses(PosePair pair, double score);

        /**
         *
         * @return The two poses yielding the best alignment of the registers ligands.
         */
        PosePair getHighestScoringPair();

        /**
         *
         * @param pose The pose to check for.
         * @return True if the register contains the @p pose.
         */
        bool containsPose(const UniquePoseIdentifier& pose);

        PosePair getHighestScoringPosePairForPose(const UniquePoseIdentifier& pose);

      private:
        LigandID m_first;
        LigandID m_second;
        unsigned m_maxSize;
        std::vector<PosePairAndScore> m_register;
    };

}  // namespace coaler::multialign
