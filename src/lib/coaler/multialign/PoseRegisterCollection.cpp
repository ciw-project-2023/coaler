/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#include "PoseRegisterCollection.hpp"

namespace coaler::multialign {

    void PoseRegisterCollection::addRegister(const PoseRegisterPtr &poseRegister) {
        LigandPair const pair(poseRegister->getFirstLigandID(), poseRegister->getSecondLigandID());

        m_registers.emplace(pair, poseRegister);
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    PairwisePoseRegisters PoseRegisterCollection::getAllRegistersForPose(const UniquePoseID &pose) const noexcept {
        PairwisePoseRegisters registersContainingPose;

        for (const auto &[ligandPair, poseRegister] : m_registers) {
            if (poseRegister->containsPose(pose)) {
                registersContainingPose.emplace(ligandPair, poseRegister);
            }
        }

        return registersContainingPose;
    }

    PairwisePoseRegisters PoseRegisterCollection::getAllRegisters() const noexcept { return m_registers; }

    /*----------------------------------------------------------------------------------------------------------------*/

}  // namespace coaler::multialign
