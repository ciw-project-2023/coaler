//
// Created by chris on 12/30/23.
//

#include "AssemblyOptimizer.hpp"

#include <spdlog/spdlog.h>

#include <coaler/embedder/ConformerEmbedder.hpp>

#include "coaler/multialign/scorer/AssemblyScorer.hpp"
#include <GraphMol/SmilesParse/SmilesWrite.h>


using namespace coaler::multialign;

void updatePoseRegisters(const LigandID ligandId, const PoseID newPose, PoseRegisterCollection &registers,
                         PairwiseAlignments &scores, const LigandVector &ligands) {
    for (const Ligand &otherLigand : ligands) {
        if (otherLigand.getID() == ligandId) {
            continue;
        }
        const LigandPair pair(ligandId, otherLigand.getID());
        for (const UniquePoseID &otherPose : otherLigand.getPoses()) {
            PosePair poses({ligandId, newPose}, otherPose);
            const double score = scores.at(poses, ligands, true);
            registers.addPoseToRegister(pair, poses, score);
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------*/

class LigandAvailabilityMapping : public std::unordered_map<LigandID, bool> {
  public:
    void setAllAvailable() {
        for (auto &pair : *this) {
            pair.second = true;
        }
    }

    /*------------------------------------------------------------------------------------------------------------*/

    void init(const LigandVector &ligands) {
        for (const Ligand &ligand : ligands) {
            this->emplace(ligand.getID(), true);
        }
    }
};

/*----------------------------------------------------------------------------------------------------------------*/

struct LigandIsAvailable {
    bool operator()(std::pair<LigandID, bool> entry) { return entry.second; }
};

/*----------------------------------------------------------------------------------------------------------------*/

LigandID get_next_missing_ligand(const LigandAlignmentAssembly &assembly, const LigandAvailabilityMapping &availability,
                                 unsigned maxLigandID) {
    auto assemblyMapping = assembly.getAssemblyMapping();
    for (LigandID id = 0; id <= maxLigandID; id++) {
        if (assemblyMapping.count(id) == 0 && availability.at(id)) {
            return id;
        }
    }
    return maxLigandID + 1;
}
/*----------------------------------------------------------------------------------------------------------------*/

OptimizerState AssemblyOptimizer::optimizeAssembly(LigandAlignmentAssembly assembly, PairwiseAlignments scores,
                                                   LigandVector ligands, PoseRegisterCollection registers,
                                                   double scoreDeficitThreshold) {
    double currentAssemblyScore = AssemblyScorer::calculateAssemblyScore(assembly, scores, ligands);

    LigandAvailabilityMapping ligandAvailable;
    ligandAvailable.init(ligands);

    // assembly optimization step
    unsigned stepCount = 0;
    while (std::any_of(ligandAvailable.begin(), ligandAvailable.end(), LigandIsAvailable())
           && stepCount < Constants::OPTIMIZER_STEP_LIMIT) {
        for (const auto &ligand : ligands) {
            if (ligand.getNumPoses() != ligand.getMolecule().getNumConformers()) {
                spdlog::error("conformer count mismatch: ligand poses {} vs. mol confs {}", ligand.getNumPoses(),
                              ligand.getMolecule().getNumConformers());
                assert(false);
            } // todo rewrite std::allof
        }
        stepCount++;
        double maxScoreDeficit = -1;
        LigandID worstLigandId = 0;
        if (assembly.getMissingLigandsCount() != 0) {
            worstLigandId = get_next_missing_ligand(assembly, ligandAvailable, ligands.size() - 1);
            if (worstLigandId == ligands.size()) {
                spdlog::debug("all invalid ligands are unavailable.");
                break;
            }
        } else {
            // no missing ligands
            for (const Ligand &ligand : ligands) {
                LigandID ligandId = ligand.getID();
                if (!ligandAvailable.at(ligand.getID())) {
                    continue;
                }
                const double ligandScoreDeficit
                    = AssemblyScorer::calculateScoreDeficitForLigand(ligandId, assembly, registers, scores, ligands);
                if (maxScoreDeficit < ligandScoreDeficit) {
                    worstLigandId = ligandId;
                    maxScoreDeficit = ligandScoreDeficit;
                }
            }
        }
        spdlog::debug("worst ligand: {} has score deficit {}", worstLigandId, maxScoreDeficit);
        if (maxScoreDeficit == 0) {
            // all pairwise alignments are optimal
            // TODO can we return this assembly and be sure its the optimum?
            // maybe once finetuning isnt score diff dependant
            break;
        }

        Ligand *worstLigand = &ligands.at(worstLigandId);
        bool ligandIsMissing = (maxScoreDeficit == -1);
        bool swappedLigandPose = false;

        // try to swap conformer of worst ligand with another existing conformer
        if (!ligandIsMissing) {
            for (const UniquePoseID &pose : worstLigand->getPoses()) {
                // check if using this pose improves currentAssembly
                if (pose.getLigandInternalPoseId() == assembly.getPoseOfLigand(worstLigandId)) {
                    continue;
                }
                LigandAlignmentAssembly assemblyCopy = assembly;
                assemblyCopy.swapPoseForLigand(worstLigandId, pose.getLigandInternalPoseId());
                double const newAssemblyScore = AssemblyScorer::calculateAssemblyScore(assemblyCopy, scores, ligands);

                if (newAssemblyScore > currentAssemblyScore) {
                    spdlog::debug("swapped for existing pose.");
                    assembly = assemblyCopy;
                    currentAssemblyScore = newAssemblyScore;
                    ligandAvailable.setAllAvailable();
                    swappedLigandPose = true;
                    break;
                }
            }
        }

        // if no improving pose can be found among existing poses, generate new ones
        if (ligandIsMissing || (!swappedLigandPose && maxScoreDeficit > scoreDeficitThreshold)) {
            spdlog::debug("generating new conformer, missing ligand = {}", ligandIsMissing);
            LigandVector alignmentTargets = {ligands.begin(), ligands.end()};

            // remove the worst ligand from targets, we only want to use all other ligands as alignment target
            for (auto target = alignmentTargets.begin(); target != alignmentTargets.end(); target++) {
                if (target->getID() == worstLigandId) {
                    alignmentTargets.erase(target);
                    break;
                }
            }

            auto newConfIDs = coaler::embedder::ConformerEmbedder::generateNewPosesForAssemblyLigand(
                (RDKit::ROMol *)worstLigand->getMoleculePtr(), alignmentTargets, assembly.getAssemblyMapping());

            if (newConfIDs.empty()) {
                spdlog::warn("no confs generated. skipping ligand {}", RDKit::MolToSmiles(worstLigand->getMolecule()));
                ligandAvailable.at(worstLigandId) = false;
                continue;
            }

            PoseID bestNewPoseID = 0;
            double bestNewAssemblyScore = 0;

            // identify new pose that yields best alignment
            for (auto iter = newConfIDs.begin(); iter != newConfIDs.end(); iter++) {
                const unsigned newPoseID = *iter;
                auto assemblyCopy = assembly;
                assemblyCopy.swapPoseForLigand(worstLigandId, newPoseID);

                // evaluate assembly
                double const newAssemblyScore = AssemblyScorer::calculateAssemblyScore(assemblyCopy, scores, ligands);
                if (newAssemblyScore > bestNewAssemblyScore) {
                    bestNewAssemblyScore = newAssemblyScore;
                    bestNewPoseID = newPoseID;
                }
            }

            if (ligandIsMissing || bestNewAssemblyScore > currentAssemblyScore) {
                currentAssemblyScore = bestNewAssemblyScore;
                spdlog::debug("ligand {} now has conformer {}.", worstLigandId, bestNewPoseID);

                // remove all (except best) new poses from ligand
                for (auto iter = newConfIDs.begin(); iter != newConfIDs.end(); iter++) {
                    if (*iter == bestNewPoseID) {
                        continue;
                    }
                    worstLigand->removePose(*iter);
                }

                // from here on we keep the new pose and adapt all containers accordingly
                updatePoseRegisters(worstLigandId, bestNewPoseID, registers, scores, ligands);
                assembly.swapPoseForLigand(worstLigandId, bestNewPoseID);
                worstLigand->addPose(bestNewPoseID);
                ligandAvailable.setAllAvailable();
                if (ligandIsMissing) {
                    assembly.decrementMissingLigandsCount();
                }
            } else {
                spdlog::debug("discarded pose. assembly score: {}", currentAssemblyScore);

                // remove all new poses from ligand
                for (auto iter = newConfIDs.begin(); iter != newConfIDs.end(); iter++) {
                    worstLigand->removePose(*iter);
                }
            }
        }
        ligands.at(worstLigandId) = *worstLigand;
        assert(ligands.at(worstLigandId).getMoleculePtr()->getNumConformers()
               == ligands.at(worstLigandId).getNumPoses());
        // set this to false in order to not immediately change this ligand again
        ligandAvailable.at(worstLigandId) = false;
    }
    spdlog::info("optimization took {} steps.", stepCount);
    OptimizerState result{currentAssemblyScore, assembly, scores, ligands, registers};
    return result;
}

/*----------------------------------------------------------------------------------------------------------------*/

OptimizerState AssemblyOptimizer::optimizeAssembly(OptimizerState &state, double scoreDeficitThreshold) {
    return optimizeAssembly(state.assembly, state.scores, state.ligands, state.registers, scoreDeficitThreshold);
}