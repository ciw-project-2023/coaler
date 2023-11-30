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
            if (registers.count(LigandPair(ligandId, otherLigand.getID())) != 0) {
                assembly.incrementMissingLigandsCount();
            }
            PosePair pair = registers.at({ligandId, otherLigand.getID()})->getHighestScoringPair();
            assert(pair.getFirst() == pose || pair.getSecond() == pose);

            UniquePoseIdentifier otherPose = pair.getFirst() == pose ? pair.getSecond() : pair.getFirst();
            assembly.insertLigandPose(otherPose.getLigandId(), otherPose.getLigandInternalPoseId());
        }
        return assembly;
    }
}  // namespace coaler::multialign