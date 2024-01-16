#pragma once

#include "PoseRegisterCollection.hpp"
#include "models/Forward.hpp"

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
