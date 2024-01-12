//
// Created by chris on 12/30/23.
//
#pragma once

#include "models/Forward.hpp"
#include "PoseRegisterCollection.hpp"

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