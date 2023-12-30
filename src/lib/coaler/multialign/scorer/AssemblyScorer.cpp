/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#include "AssemblyScorer.hpp"

#include <spdlog/spdlog.h>

#include "coaler/multialign/models/PairwiseAlignments.hpp"

namespace coaler::multialign {
    double AssemblyScorer::calculateAssemblyScore(const LigandAlignmentAssembly& assembly, PairwiseAlignments& scores,
                                                  const LigandVector& ligands) {
        double assemblyScore = 0.0;
        for (const Ligand& firstLigand : ligands) {
            for (const Ligand& secondLigand : ligands) {
                if (firstLigand.getID() >= secondLigand.getID()) {
                    continue;
                }

                PoseID const firstLigandPoseID = assembly.getPoseOfLigand(firstLigand.getID());
                PoseID const secondLigandPoseID = assembly.getPoseOfLigand(secondLigand.getID());

                // check whether assembly didnt contain one of the ligands
                if (firstLigandPoseID == std::numeric_limits<PoseID>::max()
                    || secondLigandPoseID == std::numeric_limits<PoseID>::max()) {
                    continue;
                }

                assemblyScore += getScoreInAssembly(firstLigand.getID(), secondLigand.getID(), firstLigandPoseID,
                                                    secondLigandPoseID, scores, ligands);
            }
        }
        // spdlog::info(assemblyScore);
        return assemblyScore;
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    double AssemblyScorer::calculateScoreDeficitForLigand(const LigandID ligandId, const LigandID maxLigandId,
                                                          const LigandAlignmentAssembly& assembly,
                                                          const PoseRegisterCollection& registers,
                                                          PairwiseAlignments& scores, const LigandVector& ligands) {
        PairwisePoseRegisters poseRegisters = registers.getAllRegisters();
        double scoreDeficit = 0.0;

        for (LigandID id = 0; id <= maxLigandId; id++) {
            if (id == ligandId) {
                continue;
            }
            double const scoreInAssembly = getScoreInAssembly(ligandId, id, assembly.getPoseOfLigand(ligandId),
                                                              assembly.getPoseOfLigand(id), scores, ligands);
            double const optimalScore = scores.at(poseRegisters.at(LigandPair(id, ligandId)).getHighestScoringPair());
            // TODO change this, its only a heuristic
            if (optimalScore < scoreInAssembly) {
                continue;
            }
            scoreDeficit += std::abs(optimalScore - scoreInAssembly);
        }
        return scoreDeficit;
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    double AssemblyScorer::getScoreInAssembly(LigandID firstLigandID, LigandID secondLigandID, PoseID firstPoseID,
                                              PoseID secondPoseID, PairwiseAlignments& scores,
                                              const LigandVector& ligands) {
        UniquePoseID const firstLigandPose{firstLigandID, firstPoseID};
        UniquePoseID const secondLigandPose{secondLigandID, secondPoseID};
        const PosePair pair{firstLigandPose, secondLigandPose};
        return scores.at(pair, ligands);
    }

}  // namespace coaler::multialign
