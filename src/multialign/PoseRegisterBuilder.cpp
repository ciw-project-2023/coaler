//
// Created by chris on 11/9/23.
//

#include "PoseRegisterBuilder.hpp"
#include "Ligand.hpp"
#include "PoseRegister.hpp"


namespace
{
    unsigned calculateRegisterSizeForLigand(MultiAlign::Ligand firstLigand, MultiAlign::Ligand secondLigand)
    {
        return (int) sqrt(firstLigand.getHeavyAtomSize() + secondLigand.getHeavyAtomSize());
    }
}

namespace MultiAlign {

    PairwisePoseRegisters PoseRegisterBuilder::buildPoseRegisters(
            const PairwiseAlignment &alignmentScores,
            const std::vector<Ligand>& ligands) noexcept {

        PairwisePoseRegisters registerMatrix;

        for(LigandID firstLigand = 0; firstLigand < ligands.size(); firstLigand++)
        {
            for(LigandID secondLigand = 0; secondLigand < firstLigand; secondLigand++)
            {
                if(firstLigand == secondLigand)
                {
                    continue;
                }
                unsigned size = calculateRegisterSizeForLigand(ligands.at(firstLigand), ligands.at(secondLigand));
                LigandPair currentLigandPair(firstLigand, secondLigand);
                registerMatrix.emplace(
                        currentLigandPair,std::make_shared<PoseRegister>(PoseRegister(firstLigand, secondLigand, size)));

                for(const PoseID firstLigandPose : ligands.at(firstLigand).getPoses())
                {
                    for(const PoseID secondLigandPose : ligands.at(secondLigand).getPoses())
                    {
                        registerMatrix.at(currentLigandPair)->addPoses(
                                PosePair(firstLigandPose, secondLigandPose),
                                alignmentScores.getValue(firstLigandPose, secondLigandPose));
                    }
                }
            }
        }
    }
}