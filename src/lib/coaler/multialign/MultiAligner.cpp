/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#include "MultiAligner.hpp"

#include <GraphMol/ShapeHelpers/ShapeUtils.h>
#include <omp.h>
#include <spdlog/spdlog.h>

#include <queue>
#include <utility>

#include "../embedder/ConformerEmbedder.hpp"
#include "AssemblyScorer.hpp"
#include "LigandAlignmentAssembly.hpp"
#include "StartingAssemblyGenerator.hpp"

namespace coaler::multialign {
    using AssemblyWithScore = std::pair<LigandAlignmentAssembly, double>;

    class LigandAvailabilityMapping : public std::unordered_map<LigandID, bool> {
      public:
        void setAllAvailable() {
            for (auto &pair : *this) {
                pair.second = true;
            }
        }

        /*----------------------------------------------------------------------------------------------------------------*/

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

    struct AssemblyWithScoreLess {
        bool operator()(const AssemblyWithScore &lhs, const AssemblyWithScore &rhs) {
            if (lhs.first.getMissingLigandsCount() != lhs.first.getMissingLigandsCount()) {
                return lhs.first.getMissingLigandsCount() < lhs.first.getMissingLigandsCount();
            }

            return lhs.second < rhs.second;
        }
    };

    /*----------------------------------------------------------------------------------------------------------------*/

    struct AssemblyWithScoreGreater {
        bool operator()(const AssemblyWithScore &lhs, const AssemblyWithScore &rhs) {
            if (lhs.first.getMissingLigandsCount() != lhs.first.getMissingLigandsCount()) {
                return lhs.first.getMissingLigandsCount() > lhs.first.getMissingLigandsCount();
            }

            return lhs.second > rhs.second;
        }
    };

    /*----------------------------------------------------------------------------------------------------------------*/

    MultiAligner::MultiAligner(RDKit::MOL_SPTR_VECT molecules, unsigned maxStartingAssemblies, unsigned nofThreads)

        : m_maxStartingAssemblies(maxStartingAssemblies), m_nofThreads(nofThreads) {
        assert(m_maxStartingAssemblies > 0);
        for (LigandID id = 0; id < molecules.size(); id++) {
            UniquePoseSet poses;

            for (PoseID poseId = 0; poseId < molecules.at(id)->getNumConformers(); poseId++) {
                poses.emplace(id, poseId);
            }

            m_ligands.emplace_back(*molecules.at(id), poses, id);
        }
        omp_set_num_threads(m_nofThreads);  // this sets the number of threads used for ALL subsequent parallel regions.
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    PairwiseAlignment MultiAligner::calculateAlignmentScores(const LigandVector &ligands) {
        PairwiseAlignment scores;

        // calculate number of combinations. Each pair of ligands A,B has
        // A.getNumPoses() * B.getNumPoses() many embeddings
        unsigned n = ligands.size();
        unsigned combinations = 0;
        for (unsigned id_A = 0; id_A < n; id_A++) {
            for (unsigned id_B = id_A + 1; id_B < n; id_B++) {
                combinations += ligands.at(id_A).getNumPoses() * ligands.at(id_B).getNumPoses();
            }
        }

        spdlog::info("Calculating {} combinations. This may take some time", combinations);

        for (LigandID firstMolId = 0; firstMolId < ligands.size(); firstMolId++) {
            spdlog::info("calculated {} combinations so far.", scores.size());
            for (LigandID secondMolId = firstMolId + 1; secondMolId < ligands.size(); secondMolId++) {
                unsigned nofPosesFirst = ligands.at(firstMolId).getNumPoses();
                unsigned nofPosesSecond = ligands.at(secondMolId).getNumPoses();
                omp_lock_t maplock;
                omp_init_lock(&maplock);

#pragma omp parallel for shared(maplock, ligands, scores, nofPosesFirst, nofPosesSecond, firstMolId, \
                                    secondMolId) default(none)
                for (unsigned firstMolPoseId = 0; firstMolPoseId < nofPosesFirst; firstMolPoseId++) {
                    for (unsigned secondMolPoseId = 0; secondMolPoseId < nofPosesSecond; secondMolPoseId++) {
                        const double score = getScore(ligands.at(firstMolId), ligands.at(secondMolId), firstMolPoseId,
                                                      secondMolPoseId);
                        const UniquePoseID firstPose(firstMolId, firstMolPoseId);
                        const UniquePoseID secondPose(secondMolId, secondMolPoseId);
                        omp_set_lock(&maplock);
                        scores.emplace(PosePair(firstPose, secondPose), score);
                        omp_unset_lock(&maplock);
                    }
                }
            }
        }
        spdlog::info("finished calculating pairwise alignments");
        return scores;
    }

    /*----------------------------------------------------------------------------------------------------------------*/
    MultiAlignerResult MultiAligner::alignMolecules() {
        // calculate pairwise alignments
        m_pairwiseAlignments = this->calculateAlignmentScores(m_ligands);

        spdlog::info("Mols: {} | Confs/Mol: {} | total pairwise scores: {}", m_ligands.size(),
                     m_ligands.begin()->getNumPoses(), m_pairwiseAlignments.size());
        // build pose registers
        spdlog::info("Start building pose registers.");
        m_poseRegisters = PoseRegisterBuilder::buildPoseRegisters(m_pairwiseAlignments, m_ligands, m_nofThreads);
        spdlog::info("Finish building pose registers.");

        // build starting ensembles from registers
        // AssemblyCollection assemblies;
        std::priority_queue<AssemblyWithScore, std::vector<AssemblyWithScore>, AssemblyWithScoreGreater> assemblies;
        for (const Ligand &ligand : m_ligands) {
            for (const UniquePoseID &pose : ligand.getPoses()) {
                const LigandAlignmentAssembly assembly
                    = StartingAssemblyGenerator::generateStartingAssembly(pose, m_poseRegisters, m_ligands);

                double const score = AssemblyScorer::calculateAssemblyScore(assembly, m_pairwiseAlignments, m_ligands);
                const AssemblyWithScore newAssembly = std::make_pair(assembly, score);

                // insert if queue not full or new assembly is larger that worst assembly in queue

                // TODO ensure this is called correctly
                if (assemblies.size() < m_maxStartingAssemblies) {
                    assemblies.push(newAssembly);
                    continue;
                }
                const AssemblyWithScore topAssembly = assemblies.top();
                if (AssemblyWithScoreGreater()(newAssembly, topAssembly)) {
                    assemblies.pop();
                    assemblies.push(newAssembly);
                    continue;
                }
            }
        }

        // top #m_maxStartingAssemblies are now in queue. find best scoring assembly by optimizing all
        // TODO SYMMETRY?
        // optimize all starting assemblies.
        LigandAlignmentAssembly currentBestAssembly = assemblies.top().first;  // TODO default constructor for assembly?
        double currentBestAssemblyScore
            = AssemblyScorer::calculateAssemblyScore(currentBestAssembly, m_pairwiseAlignments, m_ligands);

        // write queue content to vector to allow for parallel for
        std::vector<AssemblyWithScore> assembliesList;
        while (!assemblies.empty()) {
            assembliesList.push_back(assemblies.top());
            assemblies.pop();
        }
        spdlog::info("start optimization of {} alignment assemblies.", assembliesList.size());

        unsigned skippedAssembliesCount = 0;

        // locks for shared variables
        omp_lock_t bestAssemblyScoreLock;
        omp_init_lock(&bestAssemblyScoreLock);
        omp_lock_t bestAssemblyLock;
        omp_init_lock(&bestAssemblyLock);
        omp_lock_t skippedAssembliesCountLock;
        omp_init_lock(&skippedAssembliesCountLock);

        //#pragma omp parallel for shared(bestAssemblyScoreLock, bestAssemblyLock, skippedAssembliesCountLock, \
//                                    currentBestAssembly, currentBestAssemblyScore, assembliesList,   \
//                                    skippedAssembliesCount) default(none)

        for (unsigned assemblyID = 0; assemblyID < assembliesList.size(); assemblyID++) {
            auto [currentAssembly, currentAssemblyScore] = assembliesList.at(assemblyID);
            spdlog::info("Score before opt: {}", currentAssemblyScore);
            if (currentAssembly.getMissingLigandsCount() != 0) {
                spdlog::debug("Skip assembly because its missing ligands.");
                omp_set_lock(&skippedAssembliesCountLock);
                skippedAssembliesCount++;
                omp_unset_lock(&skippedAssembliesCountLock);
                continue;
            }

            // create a copy and uptdate it for new poses
            PairwiseAlignment pairwiseAlignments = m_pairwiseAlignments;

            LigandAvailabilityMapping ligandAvailable;
            ligandAvailable.init(m_ligands);

            // assembly optimization step
            while (std::any_of(ligandAvailable.begin(), ligandAvailable.end(), LigandIsAvailable())) {
                // determine ligand with highest score deficit TODO move to own func
                double maxScoreDeficit = 0;
                Ligand worstLigand = *m_ligands.begin();  // dummy init --> better idea?
                for (const Ligand &ligand : m_ligands) {
                    if (!ligandAvailable.at(ligand.getID())) {
                        continue;
                    }
                    const double ligandScoreDeficit = AssemblyScorer::calculateScoreDeficitForLigand(
                        ligand.getID(), m_ligands.size() - 1, currentAssembly, m_poseRegisters, pairwiseAlignments, m_ligands);
                    if (maxScoreDeficit < ligandScoreDeficit) {
                        worstLigand = ligand;
                        maxScoreDeficit = ligandScoreDeficit;
                    }
                }
                spdlog::info("max scorediff: {}", maxScoreDeficit);
                if (maxScoreDeficit == 0) {
                    // all pairwise alignments are optimal
                    // TODO can we return this assembly and be sure its the optimum?
                    break;
                }

                bool swappedLigandPose = false;
                for (const UniquePoseID &pose : worstLigand.getPoses()) {
                    // check if using this pose improves currentAssembly
                    if (pose.getLigandInternalPoseId() == currentAssembly.getPoseOfLigand(worstLigand.getID())) {
                        continue;
                    }
                    LigandAlignmentAssembly assemblyCopy = currentAssembly;
                    assemblyCopy.swapPoseForLigand(worstLigand.getID(), pose.getLigandInternalPoseId());
                    // avoid identity swap
                    double const newAssemblyScore
                        = AssemblyScorer::calculateAssemblyScore(assemblyCopy, pairwiseAlignments, m_ligands);

                    if (newAssemblyScore > currentAssemblyScore) {
                        currentAssembly = assemblyCopy;
                        currentAssemblyScore = newAssemblyScore;
                        ligandAvailable.setAllAvailable();
                        swappedLigandPose = true;
                        break;
                    }
                }

                if (!swappedLigandPose && maxScoreDeficit > 0.1) {
                    spdlog::info("generating new conformer");
                    LigandVector alignmentTargets = {m_ligands.begin(), m_ligands.end()};
                    // remove the worst ligand from targets, we only want to use all other ligands as target
                    for (auto target = alignmentTargets.begin(); target != alignmentTargets.end(); target++) {
                        if (target->getID() == worstLigand.getID()) {
                            alignmentTargets.erase(target);
                            break;
                        }
                    }

                    auto newConfIDs = coaler::embedder::ConformerEmbedder::generateNewPosesForAssemblyLigand(
                        worstLigand.getMoleculePtr().get(), alignmentTargets, currentAssembly.getAssemblyMapping());

                    if (newConfIDs.empty()) {
                        ligandAvailable.at(worstLigand.getID()) = false;
                        continue;
                    }
                    PoseID bestNewPoseID = 0;
                    double bestNewAssemblyScore = 0;

                    for (auto iter = newConfIDs.begin(); iter != newConfIDs.end(); iter++) {
                        const unsigned newPoseID = *iter;
                        auto assemblyCopy = currentAssembly;
                        assemblyCopy.swapPoseForLigand(worstLigand.getID(), newPoseID);

                        // insert pairwise scores for new pose

                        // TODO maybe use func for this? less eff but more readable
                        /*
                        for (const auto &ids : currentAssembly.getAssemblyMapping()) {
                            const LigandID otherLigandID = ids.first;
                            const PoseID otherLigandPoseID = ids.second;
                            const double score
                                = getScore(worstLigand, m_ligands.at(otherLigandID), newPoseID, otherLigandPoseID);
                            pairwiseAlignmentsCopy.emplace(
                                PosePair{{worstLigand.getID(), newPoseID}, {otherLigandID, otherLigandPoseID}}, score);
                        }*/
                        // evaluate assembly
                        double const newAssemblyScore
                            = AssemblyScorer::calculateAssemblyScore(currentAssembly, pairwiseAlignments, m_ligands);
                        if (newAssemblyScore > bestNewAssemblyScore) {
                            bestNewAssemblyScore = newAssemblyScore;
                            bestNewPoseID = newPoseID;
                        }
                    }
                    // remove all other new poses from ligand
                    for (auto iter = newConfIDs.begin(); iter != newConfIDs.end(); iter++) {
                        if (*iter == bestNewPoseID) {
                            continue;
                        }
                        m_ligands.at(worstLigand.getID()).getMolecule().removeConformer(*iter);
                    }
                    currentAssembly.swapPoseForLigand(worstLigand.getID(), bestNewPoseID);
                    // TODO add to ligand´s pose map
                    swappedLigandPose = true;
                    ligandAvailable.setAllAvailable();
                }
                if (!swappedLigandPose) {
                    ligandAvailable.at(worstLigand.getID()) = false;
                }
            }

            // check if new optimized assembly is better than current best and swap if so.
            this->ensurePairwiseAlignmentsForAssembly(currentAssembly.getAssemblyMapping());
            double const assemblyScore
                = AssemblyScorer::calculateAssemblyScore(currentAssembly, m_pairwiseAlignments, m_ligands);
            spdlog::info("Score after opt: {}", assemblyScore);
            omp_set_lock(&bestAssemblyLock);
            omp_set_lock(&bestAssemblyScoreLock);
            if (assemblyScore > currentBestAssemblyScore) {
                currentBestAssembly = currentAssembly;
                currentBestAssemblyScore = assemblyScore;
            }
            omp_unset_lock(&bestAssemblyLock);
            omp_unset_lock(&bestAssemblyScoreLock);
        }
        spdlog::info("finished alignment optimization.");
        if (skippedAssembliesCount > 0) {
            spdlog::info("Skipped a total of {} incomplete assemblies.", skippedAssembliesCount);
        }

        return {currentBestAssemblyScore, currentBestAssembly.getAssemblyMapping(), m_ligands};
    }

    [[maybe_unused]] LigandAlignmentAssembly MultiAligner::optimizeAssembly(LigandAlignmentAssembly &assembly) {}

    double MultiAligner::getScore(const Ligand &ligand1, const Ligand &ligand2, unsigned int pose1,
                                  unsigned int pose2) {
        const double distance
            = RDKit::MolShapes::tanimotoDistance(ligand1.getMolecule(), ligand2.getMolecule(), pose1, pose2);
        const double similarity = 1 - distance;
        return similarity;
    }

    void MultiAligner::ensurePairwiseAlignmentsForAssembly(const std::unordered_map<LigandID, PoseID> &assemblyIDs) {
        for (LigandID ligandId = 0; ligandId < m_ligands.size(); ligandId++) {
            const PoseID poseId = assemblyIDs.at(ligandId);
            for (LigandID otherLigandId = ligandId + 1; otherLigandId < m_ligands.size(); otherLigandId++) {
                const PoseID otherPoseId = assemblyIDs.at(ligandId);
                const UniquePoseID pose(ligandId, poseId);
                const UniquePoseID otherPose(otherLigandId, otherPoseId);
                const PosePair posePair(pose, otherPose);
                if (m_pairwiseAlignments.count(posePair) == 0) {
                    const double score
                        = getScore(m_ligands.at(ligandId), m_ligands.at(otherLigandId), poseId, otherPoseId);
                    m_pairwiseAlignments.emplace(posePair, score);
                }
            }
        }
    }

    [[maybe_unused]] void MultiAligner::addPoseToPairwiseAlignments(
        LigandID ligandId, PoseID poseId, const std::unordered_map<LigandID, PoseID> &assemblyIDs) {
        for (const auto &ids : assemblyIDs) {
            const LigandID otherLigandID = ids.first;
            if (otherLigandID == ligandId) {
                continue;
            }
            const PoseID otherPoseID = ids.second;
            const double score = getScore(m_ligands.at(ligandId), m_ligands.at(otherLigandID), poseId, otherPoseID);
            m_pairwiseAlignments.emplace(PosePair({ligandId, poseId}, {otherLigandID, otherPoseID}), score);
        }
    }

}  // namespace coaler::multialign
