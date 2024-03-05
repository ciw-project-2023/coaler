#pragma once

#include "LigandAlignmentAssembly.hpp"
#include "PoseRegisterCollection.hpp"
#include "models/Forward.hpp"

namespace coaler::multialign {
    /**
     * @brief Class to generate starting assemblies.
     */
    class StartingAssemblyGenerator {
      public:
        /**
         * @brief Generate a starting assembly.
         *
         * @param pose The pose to generate the starting assembly for.
         * @param poseCompatibilities The pose compatibilities.
         * @param ligands The ligands to generate the starting assembly for.
         * @return The generated starting assembly.
         */
        static LigandAlignmentAssembly generateStartingAssembly(UniquePoseID pose,
                                                                const PoseRegisterCollection& poseCompatibilities,
                                                                const std::vector<Ligand>& ligands);
    };
}  // namespace coaler::multialign
