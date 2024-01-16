#pragma once

#include "AssemblyOptimizer.hpp"
#include "Constants.hpp"
#include "MultiAlignerResult.hpp"
#include "coaler/core/Forward.hpp"
#include "coaler/multialign/models/Forward.hpp"

namespace coaler::multialign {

    class MultiAligner {
      public:
        explicit MultiAligner(RDKit::MOL_SPTR_VECT molecules, AssemblyOptimizer optimizer,
                              // const core::PairwiseMCSMap& pairwiseStrictMCSMap,
                              // core::PairwiseMCSMap  pairwiseRelaxedMCSMap,
                              unsigned maxStartingAssemblies = Constants::DEFAULT_NOF_STARTING_ASSEMBLIES,
                              unsigned nofThreads = Constants::DEFAULT_NOF_THREADS);

        MultiAlignerResult alignMolecules();

      private:
        static PairwiseAlignments calculateAlignmentScores(const LigandVector& ligands);

        AssemblyOptimizer m_assemblyOptimizer;

        unsigned m_maxStartingAssemblies;

        LigandVector m_ligands;
        PoseRegisterCollection m_poseRegisters;
        PairwiseAlignments m_pairwiseAlignments;

        unsigned m_threads;

        core::PairwiseMCSMap m_pairwiseStrictMcsMap;
        core::PairwiseMCSMap m_pairwiseRelaxedMcsMap;
    };

}  // namespace coaler::multialign
