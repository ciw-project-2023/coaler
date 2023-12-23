/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */
#include "PoseRegisterBuilder.hpp"

#include <spdlog/spdlog.h>

#include <iostream>

#include "PoseRegister.hpp"
#include "models/Ligand.hpp"

namespace coaler::multialign {

    PoseRegisterCollection PoseRegisterBuilder::buildPoseRegisters(const PairwiseAlignment &alignmentScores,
                                                                   const std::vector<Ligand> &ligands) noexcept {
        PairwisePoseRegisters poseRegisters;
        for (LigandID firstLigand = 0; firstLigand < ligands.size(); firstLigand++) {
            for (LigandID secondLigand = 0; secondLigand < firstLigand; secondLigand++) {
                if (firstLigand == secondLigand) {
                    continue;
                }

                unsigned size = calculateRegisterSizeForLigand(ligands.at(firstLigand), ligands.at(secondLigand));
                LigandPair currentLigandPair(firstLigand, secondLigand);

                std::shared_ptr<PoseRegister> registerPtr
                    = std::make_shared<PoseRegister>(PoseRegister(firstLigand, secondLigand, size));

                poseRegisters.emplace(currentLigandPair, registerPtr);

                for (const UniquePoseID firstLigandPose : ligands.at(firstLigand).getPoses()) {
                    for (const UniquePoseID secondLigandPose : ligands.at(secondLigand).getPoses()) {
                        PosePair pair(firstLigandPose, secondLigandPose);
                        poseRegisters.at(currentLigandPair)->addPoses(pair, alignmentScores.at(pair));
                    }
                }
            }
        }

        PoseRegisterCollection collection;
        for (const auto &reg : poseRegisters) {
            collection.addRegister(reg.second);
        }

        return collection;
    }

    unsigned PoseRegisterBuilder::calculateRegisterSizeForLigand(const Ligand &firstLigand,
                                                                 const Ligand &secondLigand) {
        //return 2* (firstLigand.getNumHeavyAtoms() + secondLigand.getNumHeavyAtoms());  // TODO find appropriate value
        double size = Constants::POSE_REGISTER_SIZE_FACTOR * firstLigand.getNumPoses() * secondLigand.getNumPoses();
        return static_cast<unsigned>(size);
    }
}  // namespace coaler::multialign
