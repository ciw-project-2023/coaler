//
// Created by chris on 11/23/23.
//

#pragma once

#include "LigandAlignmentAssembly.hpp"

namespace MultiAlign
{
    class StartingAssemblyGenerator {
        LigandAlignmentAssembly generateStartingAssembly(
                UniquePoseIdentifier Pose,
                std::shared_ptr<PoseRegister> poseCompatibilities
                );

    };
}


