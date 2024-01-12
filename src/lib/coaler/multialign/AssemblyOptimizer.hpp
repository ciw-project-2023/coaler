//
// Created by chris on 12/30/23.
//

#pragma once
#include "models/Forward.hpp"
#include "coaler/core/Forward.hpp"

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
         * @param coarseScoreThreshold Score deficits above this value will trigger the
         * generation of a new pose.
         * @param fineScoreTreshold Score deficits above this value will trigger the
         * @return The optimized assembly state
         *
         * @note You can re-call the optimization with a smaller threshold (for refinment) using the overloaded function
         */

        AssemblyOptimizer(core::PairwiseMCSMap& strictMCSMap, core::PairwiseMCSMap& relaxedMCSMap,
                          double coarseScoreThreshold, double fineScoreTreshold, int stepLimit, int threads);

        OptimizerState optimizeAssembly(LigandAlignmentAssembly assembly, PairwiseAlignments scores,
                                        LigandVector ligands, PoseRegisterCollection registers,
                                        double scoreDeficitThreshold = 0);

        /**
         * @overload
         *
         * @param state The assembly state to optimize
         * @param scoreDeficitThreshold Score deficits above this value will trigger the
         * generation of a new pose.
         * @return The optimized state.
         */
        OptimizerState fineTuneState(OptimizerState& state);

      private:
        int m_threads;
        int m_stepLimit;

        double m_coarseScoreThreshold;
        double m_fineScoreThreshold;

        core::PairwiseMCSMap& m_strictMCSMap;
        core::PairwiseMCSMap& m_relaxedMCSMap;
    };
}  // namespace coaler::multialign
