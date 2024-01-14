#include "StartingAssemblyGenerator.hpp"

namespace coaler::multialign {

    LigandAlignmentAssembly StartingAssemblyGenerator::generateStartingAssembly(
        UniquePoseID pose, const PoseRegisterCollection& poseCompatibilities, const std::vector<Ligand>& ligands) {
        LigandAlignmentAssembly assembly;

        // add already defined pose
        assembly.insertLigandPose(pose.getLigandId(), pose.getLigandInternalPoseId());
        PairwisePoseRegisters registers = poseCompatibilities.getAllRegistersForPose(pose);

        if (registers.empty()) {
            // in case that pose is in no registers
            assembly.setMissingLigandsCount(ligands.size() - 1);

            return assembly;
        }

        LigandID const ligandId = pose.getLigandId();
        for (const Ligand& otherLigand : ligands) {
            if (ligandId == otherLigand.getID()) {
                continue;
            }

            LigandPair const ligandPair(ligandId, otherLigand.getID());
            if (registers.count(ligandPair) == 0) {
                assembly.incrementMissingLigandsCount();

                continue;
            }

            PosePair const highestScoringPair = registers.at(ligandPair).getHighestScoringPosePairForPose(pose);

            auto otherPose = highestScoringPair.getFirst() == pose ? highestScoringPair.getSecond()
                                                                   : highestScoringPair.getFirst();

            assembly.insertLigandPose(otherPose.getLigandId(), otherPose.getLigandInternalPoseId());
        }

        return assembly;
    }
}  // namespace coaler::multialign
