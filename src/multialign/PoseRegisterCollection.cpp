//
// Created by chris on 11/25/23.
//

#include "PoseRegisterCollection.hpp"

namespace MultiAlign
{

    void PoseRegisterCollection::addRegister(const PoseRegisterPtr &poseRegister) {
        m_registers.emplace(LigandPair(poseRegister->getFirstLigandID(),
                                      poseRegister->getSecondLigandID()),
                           poseRegister);
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    PairwisePoseRegisters PoseRegisterCollection::getAllRegistersForPose(const UniquePoseIdentifier &pair) {
        return MultiAlign::PairwisePoseRegisters();
    }

    /*----------------------------------------------------------------------------------------------------------------*/

}