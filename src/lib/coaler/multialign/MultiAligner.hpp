/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#pragma once
#include "Forward.hpp"
#include "GraphMol/FMCS/FMCS.h"
#include "MultiAlignerResult.hpp"
#include "PoseRegister.hpp"
#include "PoseRegisterBuilder.hpp"

namespace coaler::multialign {

    class MultiAligner {
      public:
        explicit MultiAligner(RDKit::MOL_SPTR_VECT molecules,
                              unsigned maxStartingAssemblies = constants::DefaultNofStartingAssemblies,
                              unsigned nofThreads = constants::DefaultNofThreads);

        MultiAlignerResult alignMolecules();

      private:
        static PairwiseAlignment calculateAlignmentScores(const LigandVector& ligands);

        unsigned m_maxStartingAssemblies;
        std::vector<Ligand> m_ligands;
        PoseRegisterCollection m_poseRegisters;
        PairwiseAlignment m_pairwiseAlignments;
        unsigned m_nofThreads;
    };

}  // namespace coaler::multialign
