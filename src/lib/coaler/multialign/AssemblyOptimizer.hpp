//
// Created by chris on 12/30/23.
//

#pragma once
#include <coaler/multialign/models/PairwiseAlignments.hpp>
#include <coaler/core/Forward.hpp>

#include "Forward.hpp"
#include "LigandAlignmentAssembly.hpp"
#include "MultiAlignerResult.hpp"
#include "OptimizerState.hpp"
#include "PoseRegisterCollection.hpp"

namespace coaler::multialign {

    class AssemblyOptimizer {
      public:
        /**
         * Optimize an assembly using both existing conformers and
         * newly generated ones.
         *
         * @param assembly The assembly to optimize.
         * @param scores The pairwise overlap scores.
         * @param ligands The input ligands.
         * @param registers The registers for all ligand pairs.
         * @param scoreDeficitThreshold Score deficits above this value will trigger the
         * generation of a new pose.
         * @return The optimized assembly state
         *
         * @note You can re-call the optimization with a smaller threshold (for refinment) using the overloaded function
         */
        static OptimizerState optimizeAssembly(LigandAlignmentAssembly assembly, PairwiseAlignments scores,
                                               LigandVector ligands, PoseRegisterCollection registers,
                                               double scoreDeficitThreshold,
                                               const core::PairwiseMCSMap& pairwiseStrictMCSMap,
                                               const core::PairwiseMCSMap& pairwiseRelaxedMCSMap);

        /**
         * @overload
         *
         * @param state The assembly state to optimize
         * @param scoreDeficitThreshold Score deficits above this value will trigger the
         * generation of a new pose.
         * @return The optimized state.
         */
        static OptimizerState optimizeAssembly(OptimizerState& state, double scoreDeficitThreshold,
                                               const core::PairwiseMCSMap& pairwiseStrictMCSMap,
                                               const core::PairwiseMCSMap& pairwiseRelaxedMCSMap);
    };
}  // namespace coaler::multialign
