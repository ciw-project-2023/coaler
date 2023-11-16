//
// Created by chris on 11/5/23.
//

#pragma once
#include "Forward.hpp"
#include "PoseRegister.hpp"
#include "PoseRegisterBuilder.hpp"

namespace MultiAlign {

    class MultiAligner {
    public:

    private:

        std::vector<Ligand> m_ligands;
        PairwisePoseRegisters m_poseRegisters;
        PairwiseAlignment m_pairwiseAlignments;
    };

}