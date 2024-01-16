#pragma once

#include "PoseRegisterCollection.hpp"
#include "models/Forward.hpp"

namespace coaler::multialign {

    class PoseRegisterBuilder {
      public:
        // NOLINTBEGIN(readability-convert-member-functions-to-static)
        static PoseRegisterCollection buildPoseRegisters(PairwiseAlignments& alignmentScores,
                                                         const std::vector<Ligand>& ligands,
                                                         unsigned nofThreads) noexcept;

        // NOLINTEND(readability-convert-member-functions-to-static)

      private:
        static unsigned calculateRegisterSizeForLigand(const Ligand& firstLigand, const Ligand& secondLigand);
    };

}  // namespace coaler::multialign
