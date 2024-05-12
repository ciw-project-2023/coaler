#pragma once

#include "PoseRegister.hpp"
#include "models/Forward.hpp"

namespace coaler::multialign {
    /**
     * @brief Collection of PoseRegisters.
     */
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
    class PoseRegisterCollection {
      public:
        /**
         * @brief Construct a new PoseRegisterCollection object
         * @param registers The PoseRegisters to add to the collection.
         */
        void addRegister(const PoseRegister& poseRegister);

        /**
         * Get all pose registers that contain a given pose.
         * @param pose The pose to search in the register collection
         * @return A subset of all the registers; all registers containing the pose.
         */
        [[nodiscard]] PairwisePoseRegisters getAllRegistersForPose(const UniquePoseID& pose) const noexcept;

        /**
         * Get the register for a given ligand pair.
         * @param key The ligand pair to get the register for.
         * @return The register for the ligand pair.
         */
        [[nodiscard]] PoseRegisterPtr getRegisterPtr(const LigandPair& key) const noexcept;

        /**
         * Get all pose registers.
         * @return All pose registers.
         */
        [[nodiscard]] PairwisePoseRegisters getAllRegisters() const noexcept;

        /**
         * Add a pose to a register.
         * @param key The ligand pair to add the pose to.
         * @param poses The poses to add.
         * @param score The score of the alignment.
         */
        void addPoseToRegister(const LigandPair& key, const PosePair& poses, double score);

      private:
        PairwisePoseRegisters m_registers;
    };

}  // namespace coaler::multialign
