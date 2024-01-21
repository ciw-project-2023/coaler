#include "MultiAligner.hpp"

#include <omp.h>
#include <spdlog/spdlog.h>

#include <queue>
#include <utility>

#include "AssemblyOptimizer.hpp"
#include "GraphMol/FMCS/FMCS.h"
#include "LigandAlignmentAssembly.hpp"
#include "PoseRegisterBuilder.hpp"
#include "StartingAssemblyGenerator.hpp"
#include "models/AssemblyIDManager.hpp"
#include "scorer/AlignmentScorer.hpp"
#include "scorer/AssemblyScorer.hpp"

namespace coaler::multialign {
    using AssemblyWithScore = std::pair<LigandAlignmentAssembly, double>;

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
    // NOLINTBEGIN(misc-unused-parameters)
    MultiAligner::MultiAligner(RDKit::MOL_SPTR_VECT molecules, AssemblyOptimizer optimizer, core::CoreResult core,
                               unsigned maxStartingAssemblies, unsigned nofThreads)
        // NOLINTEND(misc-unused-parameters)
        : m_core(std::move(core)),
          m_maxStartingAssemblies(maxStartingAssemblies),
          m_threads(nofThreads),
          m_assemblyOptimizer(optimizer) {
        assert(m_maxStartingAssemblies > 0);

        m_ligands = LigandVector(molecules);

        // calculate pairwise alignments
        spdlog::info("start calculating pairwise alignments.");
        m_pairwiseAlignments = MultiAligner::calculateAlignmentScores(m_ligands);
        spdlog::info("finished calculating pairwise alignments.");

        // build pose registers
        spdlog::info("start building pose registers.");
        m_poseRegisters = PoseRegisterBuilder::buildPoseRegisters(m_pairwiseAlignments, m_ligands, m_threads);
        spdlog::info("finish building pose registers.");
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    PairwiseAlignments MultiAligner::calculateAlignmentScores(const LigandVector &ligands) {
        PairwiseAlignments scores;

        // calculate number of combinations. Each pair of ligands A,B has
        // A.getNumPoses() * B.getNumPoses() many embeddings
        unsigned const n = ligands.size();
        unsigned combinations = 0;
        for (unsigned idA = 0; idA < n; idA++) {
            for (unsigned idB = idA + 1; idB < n; idB++) {
                combinations += ligands.at(idA).getNumPoses() * ligands.at(idB).getNumPoses();
            }
        }

        spdlog::info("calculating {} combinations. This may take some time", combinations);

        for (LigandID firstMolId = 0; firstMolId < ligands.size(); firstMolId++) {
            spdlog::info("calculated {} combinations so far.", scores.size());
            for (LigandID secondMolId = firstMolId + 1; secondMolId < ligands.size(); secondMolId++) {
                unsigned const nofPosesFirst = ligands.at(firstMolId).getNumPoses();
                unsigned const nofPosesSecond = ligands.at(secondMolId).getNumPoses();

                omp_lock_t maplock;
                omp_init_lock(&maplock);

#pragma omp parallel for shared(maplock, ligands, scores, nofPosesFirst, nofPosesSecond, firstMolId, \
                                secondMolId) default(none)
                for (unsigned firstMolPoseId = 0; firstMolPoseId < nofPosesFirst; firstMolPoseId++) {
                    for (unsigned secondMolPoseId = 0; secondMolPoseId < nofPosesSecond; secondMolPoseId++) {
                        RDKit::RWMol const firstMol = ligands.at(firstMolId).getMolecule();
                        RDKit::RWMol const secondMol = ligands.at(secondMolId).getMolecule();
                        const double score = AlignmentScorer::calcTanimotoShapeSimilarity(
                            firstMol, secondMol, firstMolPoseId, secondMolPoseId);

                        UniquePoseID firstPose(firstMolId, firstMolPoseId);
                        UniquePoseID secondPose(secondMolId, secondMolPoseId);
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
        spdlog::info("mols: {} | confs/mol: {} | total pairwise scores: {}", m_ligands.size(),
                     m_ligands.begin()->getNumPoses(), m_pairwiseAlignments.size());

        // build starting ensembles from registers
        // AssemblyCollection assemblies;
        std::priority_queue<AssemblyWithScore, std::vector<AssemblyWithScore>, AssemblyWithScoreGreater> assemblies;

        AssemblyIDManager assemblyIdManager;
        for (const Ligand &ligand : m_ligands) {
            for (const UniquePoseID &pose : ligand.getPoses()) {
                const LigandAlignmentAssembly assembly
                    = StartingAssemblyGenerator::generateStartingAssembly(pose, m_poseRegisters, m_ligands);

                if (!assemblyIdManager.isAssemblyNew(assembly)) {
                    continue;
                }

                double const score = AssemblyScorer::calculateAssemblyScore(assembly, m_pairwiseAlignments, m_ligands);
                const AssemblyWithScore newAssembly = std::make_pair(assembly, score);

                // insert if queue not full or new assembly is larger that worst assembly in queue
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

        // write queue content to vector to allow for parallel for
        std::vector<AssemblyWithScore> assembliesList;
        while (!assemblies.empty()) {
            assembliesList.push_back(assemblies.top());
            assemblies.pop();
        }

        spdlog::info("start optimization of {} alignment assemblies.", assembliesList.size());

        unsigned skippedAssembliesCount = 0;

        // locks for shared variables
        omp_lock_t bestAssemblyLock;
        omp_init_lock(&bestAssemblyLock);
        omp_lock_t skippedAssembliesCountLock;
        omp_init_lock(&skippedAssembliesCountLock);

        OptimizerState bestAssembly{-1, {}, {}, {}, {}};

#pragma omp parallel for shared(bestAssembly, bestAssemblyLock, skippedAssembliesCount, skippedAssembliesCountLock, \
                                assembliesList) default(none)

        for (unsigned assemblyID = 0; assemblyID < assembliesList.size(); assemblyID++) {
            spdlog::debug("assembly {} has mapped Conformers for {}/{} molecules.", assemblyID,
                          assembliesList.at(assemblyID).first.getAssemblyMapping().size(), m_ligands.size());

            OptimizerState optimizedAssembly = m_assemblyOptimizer.optimizeAssembly(
                assembliesList.at(assemblyID).first, m_pairwiseAlignments, m_ligands, m_poseRegisters);

            spdlog::info("optimized assembly {}, score: {}", assemblyID, optimizedAssembly.score);

            if (optimizedAssembly.score == -1) {
                omp_set_lock(&skippedAssembliesCountLock);
                skippedAssembliesCount++;
                omp_unset_lock(&skippedAssembliesCountLock);
                continue;
            }

            omp_set_lock(&bestAssemblyLock);
            if (bestAssembly.score < optimizedAssembly.score) {
                bestAssembly = optimizedAssembly;
            }
            omp_unset_lock(&bestAssemblyLock);
        }

        // fine-tuning
        spdlog::info("fine-tuning best assembly. Score before: {}", bestAssembly.score);

        bestAssembly = m_assemblyOptimizer.fineTuneState(bestAssembly, m_core);

        spdlog::info("finished alignment optimization. Final alignment has a score of {}.", bestAssembly.score);

        if (skippedAssembliesCount > 0) {
            spdlog::info("skipped a total of {} incomplete assemblies.", skippedAssembliesCount);
        }

        MultiAlignerResult result(bestAssembly.score, bestAssembly.assembly.getAssemblyMapping(), bestAssembly.ligands);

        return result;
    }

}  // namespace coaler::multialign
