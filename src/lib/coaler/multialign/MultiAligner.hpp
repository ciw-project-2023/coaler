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
        MultiAligner(RDKit::MOL_SPTR_VECT molecules, RDKit::ROMOL_SPTR core, unsigned maxStartingAssemblies = 250);

        MultiAlignerResult alignMolecules();

      private:
        PairwiseAlignment calculateAlignmentScores(const LigandVector& ligands);

        unsigned m_maxStartingAssemblies;
        RDKit::ROMOL_SPTR m_core;
        std::vector<Ligand> m_ligands;
        PoseRegisterCollection m_poseRegisters;
        PairwiseAlignment m_pairwiseAlignments;
    };

}  // namespace coaler::multialign
