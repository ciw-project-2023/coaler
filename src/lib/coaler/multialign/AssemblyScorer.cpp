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

                const PoseID firstLigandPoseID = assembly.getPoseOfLigand(firstLigand.getID());
                const PoseID secondLigandPoseID = assembly.getPoseOfLigand(secondLigand.getID());

                // check whether assembly didnt contain one of the ligands
                if (firstLigandPoseID == std::numeric_limits<PoseID>::max()
                    || secondLigandPoseID == std::numeric_limits<PoseID>::max()) {
                    continue;
                }

                const UniquePoseID firstLigandPose{firstLigand.getID(), firstLigandPoseID};
                const UniquePoseID secondLigandPose{secondLigand.getID(), secondLigandPoseID};
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
            const UniquePoseID ligandPose(ligandId, assembly.getPoseOfLigand(ligandId));
            const UniquePoseID otherPose(id, assembly.getPoseOfLigand(id));
            const double scoreInAssembly = scores.at(PosePair(ligandPose, otherPose));
            const double optimalScore = scores.at(poseRegisters.at(LigandPair(id, ligandId))->getHighestScoringPair());

            scoreDeficit += std::abs(optimalScore - scoreInAssembly);
        }
        return scoreDeficit;
    }

}  // namespace coaler::multialign
