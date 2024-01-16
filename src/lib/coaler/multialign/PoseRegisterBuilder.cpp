#include "PoseRegisterBuilder.hpp"

#include <omp.h>
#include <spdlog/spdlog.h>

#include "Constants.hpp"
#include "PoseRegister.hpp"
#include "models/Ligand.hpp"
#include "models/PairwiseAlignments.hpp"

namespace coaler::multialign {

    // NOLINTBEGIN(misc-unused-parameters, readability-convert-member-functions-to-static)
    PoseRegisterCollection PoseRegisterBuilder::buildPoseRegisters(PairwiseAlignments &alignmentScores,
                                                                   const std::vector<Ligand> &ligands,
                                                                   unsigned nofThreads) noexcept {
        // NOLINTEND(misc-unused-parameters, readability-convert-member-functions-to-static)
        PairwisePoseRegisters poseRegisters;
        omp_lock_t poseRegistersLock;
        omp_init_lock(&poseRegistersLock);

#pragma omp parallel for default(none) shared(poseRegisters, ligands, alignmentScores, poseRegistersLock) \
    num_threads(nofThreads)
        for (LigandID firstLigand = 0; firstLigand < ligands.size(); firstLigand++) {
            for (LigandID secondLigand = 0; secondLigand < firstLigand; secondLigand++) {
                if (firstLigand == secondLigand) {
                    continue;
                }

                const unsigned size = calculateRegisterSizeForLigand(ligands.at(firstLigand), ligands.at(secondLigand));
                const LigandPair currentLigandPair(firstLigand, secondLigand);

                PoseRegister poseRegister(firstLigand, secondLigand, size);

                for (const UniquePoseID firstLigandPose : ligands.at(firstLigand).getPoses()) {
                    for (const UniquePoseID secondLigandPose : ligands.at(secondLigand).getPoses()) {
                        const PosePair pair(firstLigandPose, secondLigandPose);
                        const double score = alignmentScores.at(pair);
                        poseRegister.addPoses(pair, score);
                    }
                }

                omp_set_lock(&poseRegistersLock);
                poseRegisters.emplace(currentLigandPair, poseRegister);
                omp_unset_lock(&poseRegistersLock);
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
        // return 2* (firstLigand.getNumHeavyAtoms() + secondLigand.getNumHeavyAtoms());  // TODO find appropriate value
        double size = constants::POSE_REGISTER_SIZE_FACTOR * firstLigand.getNumPoses() * secondLigand.getNumPoses();
        return static_cast<unsigned>(size);
    }
}  // namespace coaler::multialign
