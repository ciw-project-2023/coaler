//
// Created by chris on 11/5/23.
//

#pragma once
#include "Forward.hpp"
#include "PoseRegister.hpp"
#include "PoseRegisterBuilder.hpp"
#include "GraphMol/FMCS/FMCS.h"
#include "MultiAlignerResult.hpp"
#include "../singlealign/SingleAligner.hpp"

namespace MultiAlign {

    class MultiAligner {
    public:

        MultiAligner(
                const std::vector<RDKit::RWMol>& molecules,
                RDKit::ROMol core,
                const coaler::SingleAligner& aligner);

        MultiAlignerResult alignMolecules();

    private:

        //singleAligner --> im constructior übergeben

        coaler::SingleAligner m_singleAligner;
        RDKit::ROMol m_core;
        std::vector<Ligand> m_ligands;
        PairwisePoseRegisters m_poseRegisters;
        PairwiseAlignment m_pairwiseAlignments;
    };

}