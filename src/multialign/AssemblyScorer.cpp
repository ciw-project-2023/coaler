//
// Created by chris on 11/27/23.
//

#include "AssemblyScorer.hpp"

#include <spdlog/spdlog.h>


namespace coaler::multialign
{

    double AssemblyScorer::calculateAssemblyScore(const LigandAlignmentAssembly& assembly,
                                                  const PairwiseAlignment& scores,
                                                  const LigandVector& ligands) {
        double assemblyScore = 0.0;
        for(const coaler::multialign::Ligand& firstLigand : ligands)
        {
            for(const coaler::multialign::Ligand& secondLigand : ligands)
            {
                if(firstLigand.getID() >= secondLigand.getID())
                {
                    continue;
                }
                coaler::multialign::PoseID firstLigandPoseID = assembly.getPoseOfLigand(firstLigand.getID());
                coaler::multialign::PoseID secondLigandPoseID = assembly.getPoseOfLigand(secondLigand.getID());

                if(firstLigandPoseID == std::numeric_limits<coaler::multialign::PoseID>::max()
                   || secondLigandPoseID == std::numeric_limits<coaler::multialign::PoseID>::max())
                {
                    spdlog::info("Encountered invalid PosePair during optimization.");
                    continue;
                }

                coaler::multialign::UniquePoseIdentifier firstLigandPose{
                    firstLigand.getID(), firstLigandPoseID};
                coaler::multialign::UniquePoseIdentifier secondLigandPose{
                    secondLigand.getID(), secondLigandPoseID};
                assemblyScore += scores.at(PosePair{firstLigandPose, secondLigandPose});
            }
        }
        return assemblyScore;
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    double AssemblyScorer::calculateScoreDeficitForLigand(
        const LigandID ligandId,
        const LigandID maxLigandId,
        const LigandAlignmentAssembly& assembly,
        const PoseRegisterCollection& registers,
        const PairwiseAlignment& scores)
    {
        PairwisePoseRegisters poseRegisters = registers.getAllRegisters();
        double scoreDeficit = 0.0;

        for(LigandID id = 0; id <= maxLigandId; id++)
        {
            if(id == ligandId) {continue;}
            double scoreInAssembly = scores.at(
                PosePair(
                    {id, assembly.getPoseOfLigand(id)},
                    {ligandId, assembly.getPoseOfLigand(ligandId)})
                );
            double optimalScore = scores.at(poseRegisters.at(LigandPair(id, ligandId))->getHighestScoringPair());
            scoreDeficit += (scoreInAssembly - optimalScore);
        }
        return scoreDeficit;
    }

}