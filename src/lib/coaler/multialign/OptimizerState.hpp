//
// Created by chris on 12/30/23.
//
#pragma once

#include "PoseRegisterCollection.hpp"

namespace coaler::multialign{
    struct OptimizerState{
        double score;
        LigandAlignmentAssembly assembly;
        PairwiseAlignments scores;
        LigandVector ligands;
        PoseRegisterCollection registers;

        OptimizerState& operator= (const OptimizerState& s) = default;

    };
}