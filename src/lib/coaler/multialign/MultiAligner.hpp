#pragma once

#include "AssemblyOptimizer.hpp"
#include "Constants.hpp"
#include "MultiAlignerResult.hpp"
#include "coaler/core/Forward.hpp"
#include "coaler/multialign/models/Forward.hpp"
/*!
 * @file
 * @brief Contains the MultiAligner class
 */
namespace coaler::multialign {

    /**
     * @brief The MultiAligner class is responsible for aligning multiple ligands
     * to each other.
     */
    class MultiAligner {
      public:
        /**
         * @brief Construct a new MultiAligner object
         *
         * @param molecules The molecules to align
         * @param optimizer The assembly optimizer to use
         * @param core The core result
         * @param maxStartingAssemblies The maximum number of starting assemblies to generate
         * @param nofThreads The number of threads to use
         */
        explicit MultiAligner(RDKit::MOL_SPTR_VECT molecules, AssemblyOptimizer optimizer, core::CoreResult core,
                              unsigned maxStartingAssemblies = constants::DEFAULT_NOF_STARTING_ASSEMBLIES,
                              unsigned nofThreads = constants::DEFAULT_NOF_THREADS);

        MultiAlignerResult alignMolecules();

      private:
        /**
         * @brief Generate all possible starting assemblies
         *
         * @return The generated starting assemblies
         */
        static PairwiseAlignments calculateAlignmentScores(const LigandVector& ligands);

        AssemblyOptimizer m_assemblyOptimizer;

        unsigned m_maxStartingAssemblies;

        LigandVector m_ligands;
        PoseRegisterCollection m_poseRegisters;
        PairwiseAlignments m_pairwiseAlignments;

        unsigned m_threads;

        core::CoreResult m_core;
        core::PairwiseMCSMap m_pairwiseStrictMcsMap;
        core::PairwiseMCSMap m_pairwiseRelaxedMcsMap;
    };

}  // namespace coaler::multialign
