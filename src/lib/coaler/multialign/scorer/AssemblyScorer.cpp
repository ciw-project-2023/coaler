#include "AssemblyScorer.hpp"

#include <spdlog/spdlog.h>

namespace coaler::multialign {
    double AssemblyScorer::calculateAssemblyScore(const LigandAlignmentAssembly& assembly, PairwiseAlignments& scores,
                                                  const LigandVector& ligands) {
        if (assembly.getMissingLigandsCount() > 0) {
            return 0;
        }
        double assemblyScore = 0.0;
        unsigned paircount = 0;
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
                UniquePoseID first(firstLigand.getID(), firstLigandPoseID);
                UniquePoseID second(secondLigand.getID(), secondLigandPoseID);
                assemblyScore += scores.at(PosePair{first, second}, ligands);
                paircount++;
            }
        }
        return (assemblyScore / paircount);
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    double AssemblyScorer::calculateScoreDeficitForLigand(const LigandID ligandId,
                                                          const LigandAlignmentAssembly& assembly,
                                                          const PoseRegisterCollection& registers,
                                                          PairwiseAlignments& scores, const LigandVector& ligands) {
        PairwisePoseRegisters poseRegisters = registers.getAllRegisters();
        double scoreDeficit = 0.0;

        for (LigandID otherLigandId = 0; otherLigandId < ligands.size(); otherLigandId++) {
            if (otherLigandId == ligandId) {
                continue;
            }
            PoseID firstLigandPoseId = assembly.getPoseOfLigand(ligandId);
            PoseID secondLigandPoseId = assembly.getPoseOfLigand(otherLigandId);
            UniquePoseID first(ligandId, firstLigandPoseId);
            UniquePoseID second(otherLigandId, secondLigandPoseId);
            PosePair ligandPoses = PosePair{first, second};
            double const scoreInAssembly = scores.at(ligandPoses, ligands);

            LigandPair ligandPair(ligandId, otherLigandId);
            double const optimalScore = poseRegisters.at(ligandPair).getHighestScore();

            scoreDeficit += std::abs(optimalScore - scoreInAssembly);
        }
        return scoreDeficit;
    }

    double AssemblyScorer::calculateMeanLigandDistance(LigandID ligandId, const LigandAlignmentAssembly& assembly,
                                                       PairwiseAlignments& scores, const LigandVector& ligands) {
        unsigned nofLigands = ligands.size() - 1;
        double overlapSum = 0;
        for (LigandID otherLigandId = 0; otherLigandId < ligands.size(); otherLigandId++) {
            if (otherLigandId == ligandId) {
                continue;
            }
            PoseID firstLigandPoseId = assembly.getPoseOfLigand(ligandId);
            PoseID secondLigandPoseId = assembly.getPoseOfLigand(otherLigandId);
            UniquePoseID first(ligandId, firstLigandPoseId);
            UniquePoseID second(otherLigandId, secondLigandPoseId);
            PosePair ligandPoses = PosePair{first, second};
            overlapSum += scores.at(ligandPoses, ligands);
        }
        return 1 - (overlapSum / nofLigands);
    }

}  // namespace coaler::multialign
