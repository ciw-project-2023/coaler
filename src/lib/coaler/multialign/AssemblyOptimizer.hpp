//
// Created by chris on 12/30/23.
//

#pragma once
#include <coaler/multialign/models/PairwiseAlignments.hpp>

#include "Forward.hpp"
#include "LigandAlignmentAssembly.hpp"
#include "MultiAlignerResult.hpp"
#include "OptimizerState.hpp"
#include "PoseRegisterCollection.hpp"

namespace coaler::multialign {

    class AssemblyOptimizer {
      public:
        static OptimizerState optimizeAssembly(LigandAlignmentAssembly assembly, PairwiseAlignments scores,
                                               LigandVector ligands, const PoseRegisterCollection& registers,
                                               double scoreDeficitThreshold);

        static OptimizerState optimizeAssembly(const OptimizerState& state, double scoreDeficitThreshold);
    };
}  // namespace coaler::multialign
