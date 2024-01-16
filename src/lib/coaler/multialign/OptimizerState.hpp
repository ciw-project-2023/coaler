#pragma once

#include "PoseRegisterCollection.hpp"
#include "models/Forward.hpp"

namespace coaler::multialign {

    /**
     * The OptimizerState struct represents an assembly along with everything required for the optimization.
     */
    struct OptimizerState {
        double score;
        LigandAlignmentAssembly assembly;
        PairwiseAlignments scores;
        LigandVector ligands;
        PoseRegisterCollection registers;

        OptimizerState& operator=(const OptimizerState& s) = default;
    };
}  // namespace coaler::multialign