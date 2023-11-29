//
// Created by chris on 11/25/23.
//

#include "PoseRegisterCollection.hpp"

namespace coaler::multialign
{

    void PoseRegisterCollection::addRegister(const PoseRegisterPtr &poseRegister) {
        m_registers.emplace(LigandPair(poseRegister->getFirstLigandID(),
                                      poseRegister->getSecondLigandID()),
                           poseRegister);
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    PairwisePoseRegisters PoseRegisterCollection::getAllRegistersForPose(const UniquePoseIdentifier &pose) const noexcept{

        PairwisePoseRegisters registersContainingPose;

        for(const auto& [ligandPair, poseRegister] : m_registers)
        {
            if(poseRegister->containsPose(pose))
            {
                registersContainingPose.emplace(ligandPair, poseRegister);
            }
        }

        return registersContainingPose;
    }

    PairwisePoseRegisters PoseRegisterCollection::getAllRegisters() const noexcept {
        return m_registers;
    }

    /*----------------------------------------------------------------------------------------------------------------*/

}