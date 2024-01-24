#include "AssemblyOptimizer.hpp"

#include <GraphMol/ForceFieldHelpers/UFF/UFF.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <spdlog/spdlog.h>

#include <coaler/io/OutputWriter.hpp>

#include "coaler/embedder/ConformerEmbedder.hpp"
#include "coaler/multialign/scorer/AssemblyScorer.hpp"

const unsigned SEED = 42;
const float FORCE_TOL = 0.0135;
const float RELATIVE_SCORE_THRESHOLD = 0.2;
const float ABSOLUTE_SCORE_THRESHOLD = 0.3;
const unsigned BRUTEFORCE_CONFS = 100;

namespace {
    RDKit::DGeomHelpers::EmbedParameters get_embed_params_for_optimizer_generation() {
        RDKit::DGeomHelpers::EmbedParameters params;
        params = RDKit::DGeomHelpers::srETKDGv3;
        params.optimizerForceTol = FORCE_TOL;
        params.randomSeed = SEED;
        params.useRandomCoords = true;
        params.numThreads = 1;
        params.clearConfs = false;
        return params;
    }
}  // namespace

using namespace coaler::multialign;

void update_pose_registers(const LigandID ligandId, const PoseID newPose, PoseRegisterCollection &registers,
                           PairwiseAlignments &scores, const LigandVector &ligands) {
    for (const Ligand &otherLigand : ligands) {
        if (otherLigand.getID() == ligandId) {
            continue;
        }

        const LigandPair pair(ligandId, otherLigand.getID());
        for (const UniquePoseID &otherPose : otherLigand.getPoses()) {
            const PosePair poses({ligandId, newPose}, otherPose);
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

    explicit LigandAvailabilityMapping(const LigandVector &ligands) {
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

AssemblyOptimizer::AssemblyOptimizer(coaler::core::PairwiseMCSMap &strictMCSMap,
                                     coaler::core::PairwiseMCSMap &relaxedMCSMap, embedder::ConformerEmbedder &embedder,
                                     double coarseScoreThreshold, double fineScoreThreshold, int stepLimit, int threads)
    : m_strictMCSMap(strictMCSMap),
      m_relaxedMCSMap(relaxedMCSMap),
      m_embedder(embedder),
      m_coarseScoreThreshold(coarseScoreThreshold),
      m_fineScoreThreshold(fineScoreThreshold),
      m_threads(threads),
      m_stepLimit(stepLimit) {}

/*----------------------------------------------------------------------------------------------------------------*/

std::pair<LigandID, double> get_worst_ligand_in_assembly(const LigandAlignmentAssembly &assembly,
                                                         const PoseRegisterCollection &registers,
                                                         PairwiseAlignments &scores, const LigandVector &ligands,
                                                         const LigandAvailabilityMapping &ligandAvailability) {
    double maxScoreDeficit = -1;
    LigandID worstLigandId = 0;

    if (assembly.getMissingLigandsCount() != 0) {
        worstLigandId = get_next_missing_ligand(assembly, ligandAvailability, ligands.size() - 1);
        if (worstLigandId == ligands.size()) {
            spdlog::debug("all missing ligands are unavailable.");

            return {std::numeric_limits<LigandID>::max(), maxScoreDeficit};
        }
    } else {
        // no missing ligands
        for (const Ligand &ligand : ligands) {
            const LigandID ligandId = ligand.getID();
            if (!ligandAvailability.at(ligand.getID())) {
                continue;
            }

            const double scoreDeficit
                = AssemblyScorer::calculateScoreDeficitForLigand(ligandId, assembly, registers, scores, ligands);

            if (maxScoreDeficit < scoreDeficit) {
                worstLigandId = ligandId;
                maxScoreDeficit = scoreDeficit;
            }
        }
    }

    return {worstLigandId, maxScoreDeficit};
}

/*----------------------------------------------------------------------------------------------------------------*/

LigandVector generate_alignment_targets(const LigandVector &ligands, const Ligand &ligandToAlign) {
    LigandVector alignmentTargets = {ligands.begin(), ligands.end()};

    // remove the worst ligand from targets, we only want to use all other ligands as alignment target
    auto targetsEnd = std::remove(alignmentTargets.begin(), alignmentTargets.end(), ligandToAlign);
    alignmentTargets.erase(targetsEnd, alignmentTargets.end());

    return alignmentTargets;
}

/*----------------------------------------------------------------------------------------------------------------*/

std::pair<PoseID, double> find_optimal_pose(const LigandID ligand, const std::vector<PoseID> &poses,
                                            const LigandAlignmentAssembly &assembly, PairwiseAlignments &scores,
                                            const LigandVector &ligands) {
    PoseID poseId = 0;
    double score = 0;

    // identify new pose that yields best alignment
    for (const auto newPoseId : poses) {
        auto assemblyCopy = assembly;
        assemblyCopy.swapPoseForLigand(ligand, newPoseId);

        // evaluate assembly
        double const newScore = AssemblyScorer::calculateAssemblyScore(assemblyCopy, scores, ligands);
        if (newScore > score) {
            score = newScore;
            poseId = newPoseId;
        }
    }

    return {poseId, score};
}

/*----------------------------------------------------------------------------------------------------------------*/

OptimizerState AssemblyOptimizer::optimizeAssembly(LigandAlignmentAssembly assembly, PairwiseAlignments scores,
                                                   LigandVector ligands, PoseRegisterCollection registers,
                                                   double scoreDeficitThreshold) {
    if (scoreDeficitThreshold == 0) {
        scoreDeficitThreshold = m_coarseScoreThreshold;
    }

    LigandAvailabilityMapping ligandAvailable(ligands);
    unsigned stepCount = 0;
    unsigned swapCount = 0;
    unsigned genAttempts = 0;
    unsigned genAttemptsSuccessful = 0;

    double assemblyScore = AssemblyScorer::calculateAssemblyScore(assembly, scores, ligands);

    // assembly optimization step
    auto start = std::chrono::high_resolution_clock::now();
    while (stepCount < m_stepLimit
           && std::any_of(ligandAvailable.begin(), ligandAvailable.end(), LigandIsAvailable())) {
        auto now = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::minutes>(now - start).count();

        if (duration > 3 && assembly.getMissingLigandsCount() == 0) {
            spdlog::info("optimizer run limit of 3 minuntes reached (no missing ligands).");
            break;
        }
        if (duration > 6) {
            spdlog::info("optimizer run hard limit of 6 minuntes reached (missing ligands ignored).");
            break;
        }

        assert(std::all_of(ligands.begin(), ligands.end(),
                           [](const Ligand &l) { return l.getNumPoses() == l.getMoleculePtr()->getNumConformers(); }));

        stepCount++;

        const auto [worstLigandId, maxScoreDeficit]
            = get_worst_ligand_in_assembly(assembly, registers, scores, ligands, ligandAvailable);
        spdlog::debug("worst ligand: {} has score deficit {}", worstLigandId, maxScoreDeficit);

        if (maxScoreDeficit == 0) {
            break;
        }

        if (worstLigandId == std::numeric_limits<LigandID>::max()) {
            spdlog::error(
                "Unable to generate feasible conformers using mcs method. Resorting to bruteforce conformer "
                "sampling");
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

                if (newAssemblyScore > assemblyScore) {
                    spdlog::debug("swapped for existing pose.");
                    assembly = assemblyCopy;
                    assemblyScore = newAssemblyScore;
                    if (newAssemblyScore * constants::LIGAND_AVAILABILITY_RESET_THRESHOLD > assemblyScore) {
                        ligandAvailable.setAllAvailable();
                        spdlog::debug("set available after swap");
                    }
                    swappedLigandPose = true;
                    swapCount++;
                    break;
                }
            }
        }

        const double meanDistance
            = AssemblyScorer::calculateMeanLigandDistance(worstLigandId, assembly, scores, ligands);
        if (ligandIsMissing || (!swappedLigandPose && meanDistance > scoreDeficitThreshold)) {
            spdlog::debug("generating new conformer, missing ligand = {}", ligandIsMissing);
            genAttempts++;

            LigandVector const alignmentTargets = generate_alignment_targets(ligands, *worstLigand);
            assert(alignmentTargets.size() == ligands.size() - 1);

            auto newConfIDs = coaler::embedder::ConformerEmbedder::generateNewPosesForAssemblyLigand(
                *worstLigand, alignmentTargets, assembly.getAssemblyMapping(), m_strictMCSMap, m_relaxedMCSMap,
                ligandIsMissing);

            if (newConfIDs.empty()) {
                spdlog::debug("no confs generated. skipping ligand {}", RDKit::MolToSmiles(worstLigand->getMolecule()));
                ligandAvailable.at(worstLigandId) = false;

                continue;
            }

            std::vector<std::pair<int, double>> result;
            RDKit::UFF::UFFOptimizeMoleculeConfs((RDKit::ROMol &)*worstLigand->getMoleculePtr(), result, 1);

            auto [bestNewPoseID, bestNewAssemblyScore]
                = find_optimal_pose(worstLigandId, newConfIDs, assembly, scores, ligands);

            if (ligandIsMissing || bestNewAssemblyScore > assemblyScore) {
                // from here on we keep the new pose and adapt all containers accordingly
                genAttemptsSuccessful++;
                if (ligandIsMissing
                    || bestNewAssemblyScore * constants::LIGAND_AVAILABILITY_RESET_THRESHOLD > assemblyScore) {
                    if (!ligandIsMissing) {
                        spdlog::debug("All ligands set available. Improve: {}", bestNewAssemblyScore - assemblyScore);
                    }
                    ligandAvailable.setAllAvailable();
                } else {
                    spdlog::debug("did not reset due to minor improve.");
                    ligandAvailable.at(worstLigandId) = false;
                }

                assemblyScore = bestNewAssemblyScore;
                spdlog::debug("ligand {} now has conformer {}.", worstLigandId, bestNewPoseID);

                // remove all (except best) new poses from ligand
                for (auto const confId : newConfIDs) {
                    if (confId == bestNewPoseID) {
                        continue;
                    }

                    worstLigand->removePose(confId);
                }

                update_pose_registers(worstLigandId, bestNewPoseID, registers, scores, ligands);
                assembly.swapPoseForLigand(worstLigandId, bestNewPoseID);
                worstLigand->addPose(bestNewPoseID);

                if (ligandIsMissing) {
                    assembly.decrementMissingLigandsCount();
                }

            } else {
                spdlog::debug("discarded pose. assembly score: {}", assemblyScore);

                // remove all new poses from ligand
                for (const auto confId : newConfIDs) {
                    worstLigand->removePose(confId);
                }
            }
        }

        ligands.at(worstLigandId) = *worstLigand;

        assert(worstLigand->getMoleculePtr()->getNumConformers() == ligands.at(worstLigandId).getNumPoses());

        // set this to false in order to not immediately change this ligand again
        ligandAvailable.at(worstLigandId) = false;
    }

    spdlog::info(
        "optimized assembly with score {}.\n"
        "\toptimization took {} steps:\n"
        "\t  swaps: {}\n"
        "\t  conformer generation attempts:{} ({} yielded new pose in assembly)\n",
        assemblyScore, stepCount, swapCount, genAttempts, genAttemptsSuccessful);

    return {assemblyScore, assembly, scores, ligands, registers};
}

/*----------------------------------------------------------------------------------------------------------------*/

void AssemblyOptimizer::fixWorstLigands(LigandAlignmentAssembly assembly, PairwiseAlignments scores,
                                        LigandVector ligands, PoseRegisterCollection registers) {
    spdlog::info("starting bruteforcing worst alignments in assembly.");
    double assemblyScore = AssemblyScorer::calculateAssemblyScore(assembly, scores, ligands);
    // calculating the alignment scores for all ligands separately and find average
    std::unordered_map<LigandID, double> ligandScores;
    for (Ligand &ligand : ligands) {
        unsigned ligandID = ligand.getID();
        double ligandScore = 0.0;

        for (const Ligand &secondLigand : ligands) {
            if (ligandID == secondLigand.getID()) {
                continue;
            }

            PoseID const firstLigandPoseID = assembly.getPoseOfLigand(ligandID);
            PoseID const secondLigandPoseID = assembly.getPoseOfLigand(secondLigand.getID());

            // check whether assembly didnt contain one of the ligands
            if (firstLigandPoseID == std::numeric_limits<PoseID>::max()
                || secondLigandPoseID == std::numeric_limits<PoseID>::max()) {
                continue;
            }
            const UniquePoseID first(ligandID, firstLigandPoseID);
            const UniquePoseID second(secondLigand.getID(), secondLigandPoseID);
            ligandScore += scores.at(PosePair{first, second}, ligands);
        }
        ligandScore /= (ligands.size() - 1);
        ligandScores.emplace(ligandID, ligandScore);
    }
    double ligandScoreMean = 0.0;
    for (auto &[ligandID, score] : ligandScores) {
        ligandScoreMean += score;
    }
    ligandScoreMean /= ligandScores.size();
    // bruteforce new conformers for all ligands which score is at least 10% below average
    for (Ligand &ligand : ligands) {
        LigandID ligandID = ligand.getID();
        const double currScore = ligandScores.at(ligandID);
        if (currScore + RELATIVE_SCORE_THRESHOLD * currScore < ligandScoreMean
            || currScore < ABSOLUTE_SCORE_THRESHOLD) {
            spdlog::info(
                "bruteforce: found ligand {} with below average alignment score. Starting bruteforce conformer "
                "generation.",
                RDKit::MolToSmiles(ligand.getMolecule()));
            LigandAlignmentAssembly assemblyCopy = assembly;

            // generate new conformers for ligand with fixed core coords
            const std::vector<multialign::PoseID> newPoseIDs
                = m_embedder.generateNewPosesForAssemblyLigand(ligand, BRUTEFORCE_CONFS);
            if (newPoseIDs.empty()) {
                spdlog::debug("bruteforce: no confs generated. skipping ligand {}",
                              RDKit::MolToSmiles(ligand.getMolecule()));
                continue;
            }

            std::vector<std::pair<int, double>> result;
            RDKit::UFF::UFFOptimizeMoleculeConfs((RDKit::ROMol &)*ligand.getMoleculePtr(), result, 1);

            auto [bestNewPoseID, bestNewAssemblyScore]
                = find_optimal_pose(ligandID, newPoseIDs, assemblyCopy, scores, ligands);

            spdlog::debug("bruteforce: best assembly score found with bruteforce: {}, current assembly score {}.",
                          bestNewAssemblyScore, assemblyScore);

            // add new pose to assembly if new bruteforce score is higher than old score
            if (bestNewAssemblyScore > assemblyScore) {
                assemblyScore = bestNewAssemblyScore;
                spdlog::info("bruteforce: ligand {} now has conformer {}.", ligandID, bestNewPoseID);
                // remove all (except best) new poses from ligand
                for (auto const confId : newPoseIDs) {
                    if (confId == bestNewPoseID) {
                        continue;
                    }
                    ligand.removePose(confId);
                }
                update_pose_registers(ligandID, bestNewPoseID, registers, scores, ligands);
                assembly.swapPoseForLigand(ligandID, bestNewPoseID);
                ligand.addPose(bestNewPoseID);

            } else {
                spdlog::info("bruteforce: no better conformer found for ligand {}.",
                             RDKit::MolToSmiles(ligand.getMolecule()));
                for (const auto confId : newPoseIDs) {
                    ligand.removePose(confId);
                }
            }
        }
    }
}
/*----------------------------------------------------------------------------------------------------------------*/

OptimizerState AssemblyOptimizer::fineTuneState(OptimizerState &state, const core::CoreResult &core) {
    OptimizerState optState
        = this->optimizeAssembly(state.assembly, state.scores, state.ligands, state.registers, m_fineScoreThreshold);
    this->fixWorstLigands(optState.assembly, optState.scores, optState.ligands, optState.registers);
    return optState;
}
