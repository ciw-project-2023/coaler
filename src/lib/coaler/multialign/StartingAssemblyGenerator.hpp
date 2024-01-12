/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#pragma once

#include "LigandAlignmentAssembly.hpp"
#include "PoseRegisterCollection.hpp"

#include "models/Forward.hpp"

namespace coaler::multialign {
    class StartingAssemblyGenerator {
      public:
        static LigandAlignmentAssembly generateStartingAssembly(UniquePoseID pose,
                                                                const PoseRegisterCollection& poseCompatibilities,
                                                                const std::vector<Ligand>& ligands);
    };
}  // namespace coaler::multialign
