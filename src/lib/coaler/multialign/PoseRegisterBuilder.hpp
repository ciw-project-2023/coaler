#pragma once

#include "PoseRegisterCollection.hpp"
#include "models/Forward.hpp"

namespace coaler::multialign {
    /**
     * @brief Class to build PoseRegisters.
     */
    class PoseRegisterBuilder {
      public:
        /**
         * @brief Build PoseRegisters for a set of ligands.
         *
         * @param alignmentScores The pairwise alignment scores of the ligands.
         * @param ligands The ligands to build PoseRegisters for.
         * @param nofThreads The number of threads to use.
         * @return The PoseRegisters for the ligands.
         */
        // NOLINTBEGIN(readability-convert-member-functions-to-static)
        static PoseRegisterCollection buildPoseRegisters(PairwiseAlignments& alignmentScores,
                                                         const std::vector<Ligand>& ligands,
                                                         unsigned nofThreads) noexcept;

        // NOLINTEND(readability-convert-member-functions-to-static)

      private:
        /**
         * @brief Calculate the size of the PoseRegister for a ligand.
         *
         * @param firstLigand The first ligand.
         * @param secondLigand The second ligand.
         * @return The size of the PoseRegister for the ligands.
         */
        static unsigned calculateRegisterSizeForLigand(const Ligand& firstLigand, const Ligand& secondLigand);
    };

}  // namespace coaler::multialign
