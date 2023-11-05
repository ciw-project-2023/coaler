//
// Created by chris on 11/9/23.
//
#pragma once
#include "Forward.hpp"
#include "BasicClasses/LigandPair.hpp"
#include "BasicClasses/Ligand.hpp"

namespace MultiAlign {

    class PoseRegisterBuilder {

    public:
        static PairwisePoseRegisters buildPoseRegisters(
                const PairwiseAlignment& alignmentScores,
                const std::vector<Ligand>& ligands) noexcept;

    };

}