/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#include "PoseRegisterCollection.hpp"

namespace coaler::multialign {

    void PoseRegisterCollection::addRegister(const PoseRegister& poseRegister) {
        LigandPair const pair(poseRegister.getFirstLigandID(), poseRegister.getSecondLigandID());

        m_registers.emplace(pair, poseRegister);
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    PairwisePoseRegisters PoseRegisterCollection::getAllRegistersForPose(const UniquePoseID& pose) const noexcept {
        PairwisePoseRegisters registersContainingPose;

        for (const auto& [ligandPair, poseRegister] : m_registers) {
            if (poseRegister.containsPose(pose)) {
                registersContainingPose.emplace(ligandPair, poseRegister);
            }
        }

        return registersContainingPose;
    }

    PairwisePoseRegisters PoseRegisterCollection::getAllRegisters() const noexcept { return m_registers; }

    /*----------------------------------------------------------------------------------------------------------------*/

    PoseRegisterPtr PoseRegisterCollection::getRegisterPtr(const LigandPair& key) const noexcept {
        return boost::make_shared<PoseRegister>(m_registers.at(key));
    }

    /*----------------------------------------------------------------------------------------------------------------*/

}  // namespace coaler::multialign
