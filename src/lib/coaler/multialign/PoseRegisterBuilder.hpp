/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */
#pragma once
#include "Forward.hpp"
#include "PoseRegisterCollection.hpp"
#include "models/Ligand.hpp"
#include "models/LigandPair.hpp"
#include "models/PairwiseAlignments.hpp"

namespace coaler::multialign {

    class PoseRegisterBuilder {
      public:
        static PoseRegisterCollection buildPoseRegisters(PairwiseAlignments& alignmentScores,
                                                         const std::vector<Ligand>& ligands,
                                                         unsigned nofThreads) noexcept;

      private:
        static unsigned calculateRegisterSizeForLigand(const Ligand& firstLigand, const Ligand& secondLigand);
    };

}  // namespace coaler::multialign
