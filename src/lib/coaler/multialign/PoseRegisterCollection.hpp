/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#pragma once
#include "Forward.hpp"
#include "PoseRegister.hpp"
#include "models/Forward.hpp"

namespace coaler::multialign {

    class PoseRegisterCollection {
      public:
        // TODO prefer init in constructor.
        /**
         *
         * @param poseRegister PoseRegister to add to the collection.
         */
        void addRegister(const PoseRegisterPtr& poseRegister);

        /**
         * Get all pose registers that contain a given pose.
         * @param pose The pose to search in the register collection
         * @return A subset of all the registers; all registers containing the pose.
         */
        PairwisePoseRegisters getAllRegistersForPose(const UniquePoseID& pose) const noexcept;

        PairwisePoseRegisters getAllRegisters() const noexcept;

      private:
        PairwisePoseRegisters m_registers;
    };

}  // namespace coaler::multialign
