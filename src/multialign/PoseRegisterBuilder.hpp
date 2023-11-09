//
// Created by chris on 11/9/23.
//
#pragma once
#include "Forward.hpp"
#include "LigandPair.hpp"
#include "Ligand.hpp"

namespace MultiAlign {

    using PairwisePoseRegisters = std::unordered_map<LigandPair, PoseRegisterPtr, LigandPairHash>;

    class PoseRegisterBuilder {

        static PairwisePoseRegisters buildPoseRegisters(const PairwiseAlignment& alignmentScores,
                                                       const std::vector<Ligand>& ligands) noexcept;

    };

}