//
// Created by chris on 11/5/23.
//

#pragma once
#include "Forward.hpp"
#include "PoseRegister.hpp"
#include "PoseRegisterBuilder.hpp"
#include "GraphMol/FMCS/FMCS.h"
#include "MultiAlignerResult.hpp"

namespace MultiAlign {

    class MultiAligner {
    public:

        MultiAligner(
                const std::vector<RDKit::RWMol>& molecules,
                const RDKit::MCSResult& core);

        MultiAlignerResult alignMolecules();

    private:

        std::vector<Ligand> m_ligands;
        PairwisePoseRegisters m_poseRegisters;
        PairwiseAlignment m_pairwiseAlignments;
    };

}