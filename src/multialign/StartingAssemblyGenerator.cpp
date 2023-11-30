//
// Created by chris on 11/23/23.
//

#include "StartingAssemblyGenerator.hpp"

namespace coaler::multialign {

    LigandAlignmentAssembly StartingAssemblyGenerator::generateStartingAssembly(
        UniquePoseIdentifier pose, const PoseRegisterCollection& poseCompatibilities,
        const std::vector<Ligand>& ligands) {
        LigandAlignmentAssembly assembly;

        // add already defined pose
        assembly.insertLigandPose(pose.getLigandId(), pose.getLigandInternalPoseId());
        PairwisePoseRegisters registers = poseCompatibilities.getAllRegistersForPose(pose);
        LigandID ligandId = pose.getLigandId();

        for (const Ligand& otherLigand : ligands) {
            if (ligandId == otherLigand.getID()){
                continue;
            }
            LigandPair ligandPair(ligandId, otherLigand.getID());
            if (registers.count(ligandPair) == 0) {
                assembly.incrementMissingLigandsCount();
            }
            PosePair highestScoringPair = registers.at(ligandPair)->getHighestScoringPosePairForPose(pose);

//            UniquePoseIdentifier otherPose =
 //               highestScoringPair.getFirst() == pose ? highestScoringPair.getSecond() : highestScoringPair.getFirst();
            UniquePoseIdentifier otherPose;
            if(highestScoringPair.getFirst() == pose)
            {
                otherPose = highestScoringPair.getSecond();
            } else {
                otherPose = highestScoringPair.getFirst();
            }

            assembly.insertLigandPose(otherPose.getLigandId(), otherPose.getLigandInternalPoseId());
        }
        return assembly;
    }
}  // namespace coaler::multialign