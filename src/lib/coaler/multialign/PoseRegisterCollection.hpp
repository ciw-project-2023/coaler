#pragma once

#include "PoseRegister.hpp"
#include "models/Forward.hpp"

namespace coaler::multialign {

    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
    class PoseRegisterCollection {
      public:
        /**
         * @param poseRegister PoseRegister to add to the collection.
         */
        void addRegister(const PoseRegister& poseRegister);

        /**
         * Get all pose registers that contain a given pose.
         * @param pose The pose to search in the register collection
         * @return A subset of all the registers; all registers containing the pose.
         */
        [[nodiscard]] PairwisePoseRegisters getAllRegistersForPose(const UniquePoseID& pose) const noexcept;

        [[nodiscard]] PoseRegisterPtr getRegisterPtr(const LigandPair& key) const noexcept;

        [[nodiscard]] PairwisePoseRegisters getAllRegisters() const noexcept;

        void addPoseToRegister(const LigandPair& key, const PosePair& poses, double score);

      private:
        PairwisePoseRegisters m_registers;
    };

}  // namespace coaler::multialign
