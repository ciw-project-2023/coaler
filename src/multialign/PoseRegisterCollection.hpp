//
// Created by chris on 11/25/23.
//

#pragma once
#include "BasicClasses/Forward.hpp"
#include "Forward.hpp"
#include "PoseRegister.hpp"

namespace MultiAlign
{

class PoseRegisterCollection {
public:

    void addRegister(const PoseRegisterPtr& poseRegister);

    PairwisePoseRegisters getAllRegistersForPose(const UniquePoseIdentifier& pose);

private:

    PairwisePoseRegisters m_registers;
};

}

