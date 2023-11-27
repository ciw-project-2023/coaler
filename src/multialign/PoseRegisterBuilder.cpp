//
// Created by chris on 11/9/23.
//

#include "PoseRegisterBuilder.hpp"
#include "BasicClasses/Ligand.hpp"
#include "PoseRegister.hpp"
#include <spdlog/spdlog.h>
#include <iostream>

namespace
{
    unsigned calculateRegisterSizeForLigand(MultiAlign::Ligand firstLigand, MultiAlign::Ligand secondLigand)
    {
        return (int) sqrt(firstLigand.getHeavyAtomSize() + secondLigand.getHeavyAtomSize());
    }
}

namespace MultiAlign {

    PoseRegisterCollection PoseRegisterBuilder::buildPoseRegisters(
            const PairwiseAlignment &alignmentScores,
            const std::vector<Ligand>& ligands) noexcept {

        PairwisePoseRegisters poseRegisters;
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
                poseRegisters.emplace(
                        currentLigandPair,std::make_shared<PoseRegister>(PoseRegister(firstLigand, secondLigand, size)));

                for(const UniquePoseIdentifier firstLigandPose : ligands.at(firstLigand).getPoses())
                {
                    for(const UniquePoseIdentifier secondLigandPose : ligands.at(secondLigand).getPoses())
                    {
                        PosePair pair(firstLigandPose, secondLigandPose);
                        poseRegisters.at(currentLigandPair)->addPoses(
                                pair,
                                alignmentScores.at(pair));
                    }
                }
            }
        }
        PoseRegisterCollection collection;
        for(const auto& reg : poseRegisters)
        {
            collection.addRegister(reg.second);
        }
        return collection;
    }
}