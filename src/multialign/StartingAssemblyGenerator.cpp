//
// Created by chris on 11/23/23.
//

#include "StartingAssemblyGenerator.hpp"

namespace MultiAlign
{

    LigandAlignmentAssembly StartingAssemblyGenerator::generateStartingAssembly(
            UniquePoseIdentifier pose,
            const PoseRegisterCollection& poseCompatibilities,
            const std::vector<Ligand>& ligands)
    {
        LigandAlignmentAssembly assembly;

        //add already defined pose
        assembly.insertLigandPose(pose.m_ligandId, pose.m_ligandInternalPoseId);
        PairwisePoseRegisters registers = poseCompatibilities.getAllRegistersForPose(pose);
        LigandID ligandId = pose.m_ligandId;

        for(const Ligand& otherLigand : ligands)
        {
            if(registers.count(LigandPair(ligandId, otherLigand.getID())) != 0)
            {
                assembly.incrementMissingLigandsCount();
            }
            PosePair pair = registers.at({ligandId, otherLigand.getID()})->getHighestScoringPair();
            assert(pair.getFirst() == pose || pair.getSecond() == pose);

            UniquePoseIdentifier otherPose = pair.getFirst() == pose? pair.getSecond() : pair.getFirst();
            assembly.insertLigandPose(otherPose.m_ligandId, otherPose.m_ligandInternalPoseId);
        }
        return assembly;
    }
}