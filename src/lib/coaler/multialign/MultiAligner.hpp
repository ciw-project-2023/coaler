/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#pragma once
#include "../singlealign/SingleAligner.hpp"
#include "Forward.hpp"
#include "GraphMol/FMCS/FMCS.h"
#include "MultiAlignerResult.hpp"
#include "PoseRegister.hpp"
#include "PoseRegisterBuilder.hpp"

namespace coaler::multialign {

    class MultiAligner {
      public:
        MultiAligner(const std::vector<RDKit::RWMol*>& molecules, RDKit::ROMol core,
                     const coaler::SingleAligner& aligner, unsigned maxStartingAssemblies = 250);

        MultiAlignerResult alignMolecules();

      private:
        PairwiseAlignment calculateAlignmentScores(const LigandVector& ligands, const RDKit::ROMol& core);

        unsigned m_maxStartingAssemblies;
        coaler::SingleAligner m_singleAligner;
        RDKit::ROMol m_core;
        std::vector<Ligand> m_ligands;
        PoseRegisterCollection m_poseRegisters;
        PairwiseAlignment m_pairwiseAlignments;
    };

}  // namespace coaler::multialign
