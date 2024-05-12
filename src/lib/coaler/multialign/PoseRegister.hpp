#pragma once

#include <unordered_map>

#include "boost/shared_ptr.hpp"
#include "models/Forward.hpp"

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
        void addPoses(PosePair pair, double score);

        /**
         * @return The two poses yielding the best alignment of the registers ligands.
         */
        [[nodiscard]] PosePair getHighestScoringPair() const noexcept;

        /**
         * @return The overlap score of the best scoring pose pair.
         */
        [[nodiscard]] double getHighestScore() const noexcept;

        /**
         *
         * @param pose The pose to search for
         * @return The highest scoring entry in the register that is composed by @p pose
         */
        PosePair getHighestScoringPosePairForPose(const UniquePoseID& pose);

        /**
         *
         * @param pose The pose to check for.
         * @return True if the register contains the @p pose.
         */
        [[nodiscard]] bool containsPose(const UniquePoseID& pose) const;

      private:
        void updateLowest();
        void removeAndUpdateLowest();
        void updateHighest(const PosePairAndScore& insertedPair);

        LigandID m_first;
        LigandID m_second;
        unsigned m_maxSize;
        std::vector<PosePairAndScore> m_register;
        PosePairAndScore m_lowest;
        PosePairAndScore m_highest;
    };

    using PoseRegisterPtr = boost::shared_ptr<PoseRegister>;
    using PairwisePoseRegisters = std::unordered_map<LigandPair, PoseRegister, LigandPairHash>;
}  // namespace coaler::multialign
