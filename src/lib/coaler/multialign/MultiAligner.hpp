/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#pragma once
#include <coaler/core/Forward.hpp>
#include <coaler/multialign/models/PairwiseAlignments.hpp>

#include "Forward.hpp"
#include "GraphMol/FMCS/FMCS.h"
#include "LigandAlignmentAssembly.hpp"
#include "MultiAlignerResult.hpp"
#include "PoseRegister.hpp"
#include "PoseRegisterBuilder.hpp"

namespace coaler::multialign {

    class MultiAligner {
      public:
        explicit MultiAligner(RDKit::MOL_SPTR_VECT molecules,
                              // const core::PairwiseMCSMap& pairwiseStrictMCSMap,
                              // core::PairwiseMCSMap  pairwiseRelaxedMCSMap,
                              unsigned maxStartingAssemblies = Constants::DEFAULT_NOF_STARTING_ASSEMBLIES,
                              unsigned nofThreads = Constants::DEFAULT_NOF_THREADS);

        MultiAlignerResult alignMolecules();

      private:
        static PairwiseAlignments calculateAlignmentScores(const LigandVector& ligands);

        unsigned m_maxStartingAssemblies;
        std::vector<Ligand> m_ligands;
        PoseRegisterCollection m_poseRegisters;
        PairwiseAlignments m_pairwiseAlignments;
        unsigned m_nofThreads;
        core::PairwiseMCSMap m_pairwiseStrictMcsMap;
        core::PairwiseMCSMap m_pairwiseRelaxedMcsMap;
    };

}  // namespace coaler::multialign
