//
// Created by chris on 11/9/23.
//
#pragma once
#include "Forward.hpp"
#include "BasicClasses/LigandPair.hpp"
#include "BasicClasses/Ligand.hpp"
#include "PoseRegisterCollection.hpp"


namespace coaler::multialign {

    class PoseRegisterBuilder {

    public:
        static PoseRegisterCollection buildPoseRegisters(
                const PairwiseAlignment& alignmentScores,
                const std::vector<Ligand>& ligands) noexcept;

    };

}