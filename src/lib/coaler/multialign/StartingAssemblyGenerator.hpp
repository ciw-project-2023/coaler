//
// Created by chris on 11/23/23.
//

#pragma once

#include "LigandAlignmentAssembly.hpp"
#include "PoseRegisterCollection.hpp"

namespace coaler::multialign {
    class StartingAssemblyGenerator {
      public:
        static LigandAlignmentAssembly generateStartingAssembly(UniquePoseID pose,
                                                                const PoseRegisterCollection& poseCompatibilities,
                                                                const std::vector<Ligand>& ligands);
    };
}  // namespace coaler::multialign
