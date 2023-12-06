/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#include "AssemblyScorer.hpp"

#include <spdlog/spdlog.h>

namespace coaler::multialign {
    double AssemblyScorer::calculateAssemblyScore(const LigandAlignmentAssembly& assembly,
                                                  const PairwiseAlignment& scores, const LigandVector& ligands) {
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

                UniquePoseID const firstLigandPose{firstLigand.getID(), firstLigandPoseID};
                UniquePoseID const secondLigandPose{secondLigand.getID(), secondLigandPoseID};
                assemblyScore += scores.at(PosePair{firstLigandPose, secondLigandPose});
            }
        }
        // spdlog::info(assemblyScore);
        return assemblyScore;
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    double AssemblyScorer::calculateScoreDeficitForLigand(const LigandID ligandId, const LigandID maxLigandId,
                                                          const LigandAlignmentAssembly& assembly,
                                                          const PoseRegisterCollection& registers,
                                                          const PairwiseAlignment& scores) {
        PairwisePoseRegisters poseRegisters = registers.getAllRegisters();
        double scoreDeficit = 0.0;

        for (LigandID id = 0; id <= maxLigandId; id++) {
            if (id == ligandId) {
                continue;
            }
            UniquePoseID ligandPose(ligandId, assembly.getPoseOfLigand(ligandId));
            UniquePoseID otherPose(id, assembly.getPoseOfLigand(id));
            double const scoreInAssembly = scores.at(PosePair(ligandPose, otherPose));
            double const optimalScore = scores.at(poseRegisters.at(LigandPair(id, ligandId))->getHighestScoringPair());

            scoreDeficit += std::abs(optimalScore - scoreInAssembly);
        }
        return scoreDeficit;
    }

}  // namespace coaler::multialign
