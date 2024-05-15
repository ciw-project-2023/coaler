#pragma once

#include "coaler/multialign/LigandAlignmentAssembly.hpp"
#include "coaler/multialign/PoseRegisterCollection.hpp"
#include "coaler/multialign/models/PairwiseAlignments.hpp"

namespace coaler::multialign {
    /**
     * @brief Class to score an assembly.
     */
    class AssemblyScorer {
      public:
        /**
         * @brief calculate the score of an assembly.
         * Scoring method: sum of all pairwise overlaps.
         * @param assembly The assembly to score.
         * @param scores The pairwise overlap scores.
         * @param ligands The ligands the assembly contains.
         * @return The score of the assembly.
         */
        static double calculateAssemblyScore(const coaler::multialign::LigandAlignmentAssembly& assembly,
                                             coaler::multialign::PairwiseAlignments& scores,
                                             const coaler::multialign::LigandVector& ligands);

        /**
         * @brief Calculate the score deficit of some Ligand conformer.
         *
         * The score defecit is a measure of how bad a conformer is in an assembly in comparison to the best conformer
         * of a ligand.
         *
         * @param ligandId ID of the ligand to asses.
         * @param assembly The assembly in whose context the ligand conformer is scored.
         * @param registers The pairwise pose registers of the ligands.
         * @param scores The pairwise overlap scores of the ligands´ conformers
         * @param ligands The input ligands.
         * @return The score deficit of the ligands conformer.
         */
        static double calculateScoreDeficitForLigand(LigandID ligandId, const LigandAlignmentAssembly& assembly,
                                                     const PoseRegisterCollection& registers,
                                                     PairwiseAlignments& scores, const LigandVector& ligands);

        /**
         * @brief Calculate the mean distance of a ligand to all other ligands in the assembly.
         * @param ligandId
         * @param assembly
         * @param scores
         * @param ligands
         * @return The mean distance of a ligand to all other ligands in the assembly.
         */
        static double calculateMeanLigandDistance(LigandID ligandId, const LigandAlignmentAssembly& assembly,
                                                  PairwiseAlignments& scores, const LigandVector& ligands);
    };

}  // namespace coaler::multialign
