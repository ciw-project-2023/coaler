#pragma once
#include "LigandAlignmentAssembly.hpp"
#include "MultiAlignerResult.hpp"
#include "OptimizerState.hpp"
#include "PoseRegisterCollection.hpp"
#include "coaler/core/Forward.hpp"
#include "coaler/embedder/Forward.hpp"
#include "models/Forward.hpp"

namespace coaler::multialign {

    class AssemblyOptimizer {
      public:
        /**
         * Optimize an assembly using both existing conformers and
         * newly generated ones.
         *
         * @param assembly The assembly to optimize
         * @param scores The pairwise overlap scores
         * @param ligands The input ligands
         * @param registers The registers for all ligand pairs
         * @param coarseScoreThreshold Score deficits above this value will trigger the
         * generation of a new pose
         * @param fineScoreTreshold Score deficits above this value will trigger the
         * @return The optimized assembly state
         *
         * @note You can re-call the optimization with a smaller threshold (for refinment) using the overloaded function
         */
        // NOLINTBEGIN(readability-inconsistent-declaration-parameter-name)
        AssemblyOptimizer(core::PairwiseMCSMap& strictMCSMap, core::PairwiseMCSMap& relaxedMCSMap,
                          embedder::ConformerEmbedder& embedder, double coarseScoreThreshold, double fineScoreTreshold,
                          int stepLimit, int threads);
        // NOLINTEND(readability-inconsistent-declaration-parameter-name)

        /**
         * @brief Optimize an assembly using both existing conformers and newly generated ones.
         *
         * @param assembly The assembly to optimize
         * @param scores The pairwise overlap scores
         * @param ligands The input ligands
         * @param registers The registers for all ligand pairs
         * @param scoreDeficitThreshold Score deficits above this value will trigger the
         * generation of a new pose
         * @return The optimized state
         */
        OptimizerState optimizeAssembly(LigandAlignmentAssembly assembly, PairwiseAlignments scores,
                                        LigandVector ligands, PoseRegisterCollection registers,
                                        double scoreDeficitThreshold = 0);

        /**
         * @overload
         *
         * @param state The assembly state to optimize
         * @param scoreDeficitThreshold Score deficits above this value will trigger the
         * generation of a new pose
         * @return The optimized state
         */
        OptimizerState fineTuneState(OptimizerState& state, const core::CoreResult& core);

      private:
        /**
         * Using fixed core conformer generation for ligands with a below average alignment score in the best assembly
         * @param assembly best assembly found
         * @param scores The pairwise overlap scores
         * @param ligands The input ligands
         * @param registers The registers for all ligand pairs
         * @param core core of all input molecules
         */
        void fixWorstLigands(LigandAlignmentAssembly assembly, PairwiseAlignments scores, LigandVector ligands,
                             PoseRegisterCollection registers);

        int m_threads;
        int m_stepLimit;

        double m_coarseScoreThreshold;
        double m_fineScoreThreshold;

        core::PairwiseMCSMap& m_strictMCSMap;
        core::PairwiseMCSMap& m_relaxedMCSMap;

        embedder::ConformerEmbedder& m_embedder;
    };
}  // namespace coaler::multialign
