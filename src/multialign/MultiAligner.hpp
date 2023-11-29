//
// Created by chris on 11/5/23.
//

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
        MultiAligner(const std::vector<RDKit::RWMol>& molecules, RDKit::ROMol core,
                     const coaler::SingleAligner& aligner, unsigned maxStartingAssemblies);

        MultiAlignerResult alignMolecules();

      private:
        unsigned m_maxStartingAssemblies;
        coaler::SingleAligner m_singleAligner;
        RDKit::ROMol m_core;
        std::vector<Ligand> m_ligands;
        PoseRegisterCollection m_poseRegisters;
        PairwiseAlignment m_pairwiseAlignments;
    };

}  // namespace coaler::multialign