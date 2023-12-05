/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#include "MultiAligner.hpp"

#include <spdlog/spdlog.h>

#include <queue>
#include <utility>

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

        void init(const LigandVector& ligands) {
            for(const Ligand& ligand : ligands)
            {
                this->emplace(ligand.getID(), true);
            }
        }
    };

    struct LigandIsAvailable {
        bool operator()(std::pair<LigandID, bool> entry) { return entry.second; }
    };

    struct AssemblyWithScoreLess {
        bool operator()(const AssemblyWithScore &lhs, const AssemblyWithScore &rhs) {
            if (lhs.first.getMissingLigandsCount() != lhs.first.getMissingLigandsCount()) {
                return lhs.first.getMissingLigandsCount() < lhs.first.getMissingLigandsCount();
            }

            return lhs.second < rhs.second;
        }
    };

    struct AssemblyWithScoreGreater {
        bool operator()(const AssemblyWithScore &lhs, const AssemblyWithScore &rhs) {
            if (lhs.first.getMissingLigandsCount() != lhs.first.getMissingLigandsCount()) {
                return lhs.first.getMissingLigandsCount() > lhs.first.getMissingLigandsCount();
            }

            return lhs.second > rhs.second;
        }
    };

    MultiAligner::MultiAligner(const std::vector<RDKit::RWMol*> &molecules, RDKit::ROMol core,
                               const coaler::SingleAligner &aligner, unsigned maxStartingAssemblies)

        : m_core(std::move(core)), m_singleAligner(aligner), m_maxStartingAssemblies(maxStartingAssemblies) {
        assert(m_maxStartingAssemblies > 0);
        for (LigandID id = 0; id < molecules.size(); id++) {
            UniquePoseSet poses;

            for (PoseID poseId = 0; poseId < molecules.at(id)->getNumConformers(); poseId++) {
                poses.emplace(id, poseId);
            }

            m_ligands.emplace_back(*molecules.at(id), poses, id);
        }
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    PairwiseAlignment MultiAligner::calculateAlignmentScores(const LigandVector &ligands, const RDKit::ROMol &core) {
        PairwiseAlignment scores;
        unsigned n = ligands.size();
        unsigned m = ligands.at(0).getNumPoses();
        unsigned combinations = (int) (0.5 * ((n-1)*(n-1)*m*m));
        unsigned calced = 0;
        for (LigandID firstMolId = 0; firstMolId < ligands.size(); firstMolId++) {
            for (LigandID secondMolId = firstMolId + 1; secondMolId < ligands.size(); secondMolId++) {
                unsigned nofPosesFirst = ligands.at(firstMolId).getNumPoses();
                unsigned nofPosesSecond = ligands.at(secondMolId).getNumPoses();

                for (unsigned firstMolPoseId = 0; firstMolPoseId < nofPosesFirst; firstMolPoseId++) {
                    for (unsigned secondMolPoseId = 0; secondMolPoseId < nofPosesSecond; secondMolPoseId++) {
                        RDKit::RWMol const firstMol = ligands.at(firstMolId).getMolecule();
                        RDKit::RWMol const secondMol = ligands.at(secondMolId).getMolecule();

                        double score = m_singleAligner.align_molecules_kabsch(firstMol, secondMol, firstMolPoseId,
                                                                              secondMolPoseId, core);
                        calced++;
                        UniquePoseID firstPose(firstMolId, firstMolPoseId);
                        UniquePoseID secondPose(secondMolId, secondMolPoseId);

                        scores.emplace(PosePair(firstPose, secondPose), score);
                    }
                }
                if(calced % 1000 == 0)
                {
                    spdlog::info("Calculated {}/{} combinations", calced, combinations);
                }
            }
        }
        spdlog::info("finished calculating pairwise alignments");
        return scores;
    }

    /*----------------------------------------------------------------------------------------------------------------*/
    MultiAlignerResult MultiAligner::alignMolecules() {
        // calculate pairwise alignments
        m_pairwiseAlignments = this->calculateAlignmentScores(m_ligands, m_core);

        //build pose registers
        m_poseRegisters = PoseRegisterBuilder::buildPoseRegisters(m_pairwiseAlignments, m_ligands);

        // build starting ensembles from registers
        // AssemblyCollection assemblies;
        std::priority_queue<AssemblyWithScore, std::vector<AssemblyWithScore>, AssemblyWithScoreGreater> assemblies;

        for (const Ligand &ligand : m_ligands) {
            for (const UniquePoseID &pose : ligand.getPoses()) {
                LigandAlignmentAssembly assembly
                    = StartingAssemblyGenerator::generateStartingAssembly(pose, m_poseRegisters, m_ligands);

                double const score = AssemblyScorer::calculateAssemblyScore(assembly, m_pairwiseAlignments, m_ligands);
                AssemblyWithScore newAssembly = std::make_pair(assembly, score);

                // insert if queue no full or new assembly is larger that worst assembly in queue

                // TODO ensure this is called correctly
                if (assemblies.size() < m_maxStartingAssemblies){
                    assemblies.push(newAssembly);
                    continue;
                }
                AssemblyWithScore topAssembly = assemblies.top();
                if(AssemblyWithScoreGreater()(newAssembly,topAssembly)) {
                    assemblies.pop();
                    assemblies.push(newAssembly);
                    continue;
                }
            }
        }

        // top #m_maxStartingAssemblies are now in queue. find best scoring assembly by optimizing all
        //TODO SYMMETRIE?
        // optimize all starting assemblies.
        LigandAlignmentAssembly currentBestAssembly = assemblies.top().first;  // TODO default constructor for assembly?
        double currentBestScore
            = AssemblyScorer::calculateAssemblyScore(currentBestAssembly, m_pairwiseAlignments, m_ligands);

        while (!assemblies.empty()) {
            LigandAlignmentAssembly currentAssembly = assemblies.top().first;
            spdlog::info("remaining assemblies to be opimized: {}", assemblies.size());
            double currentAssemblyScore = assemblies.top().second;
            spdlog::info("Score before opt: {}", currentAssemblyScore);
            assemblies.pop();
            if(currentAssembly.getMissingLigandsCount() != 0)
            {
                spdlog::info("Skip assembly because its missing ligands.");
                continue;
            }

            LigandAvailabilityMapping ligandAvailable;
            ligandAvailable.init(m_ligands);

            //TODO dont recalc score but use assemblies score

            while (std::any_of(ligandAvailable.begin(), ligandAvailable.end(), LigandIsAvailable())) {
                double maxScoreDeficit = 0;
                for(const Ligand& ligand : m_ligands)
                {
                    double ligandScoreDeficit =
                        AssemblyScorer::calculateScoreDeficitForLigand(ligand, currentAssembly);
                }
                for (const Ligand &ligand : m_ligands) {
                    if (!ligandAvailable.at(ligand.getID())) {
                        continue;
                    }
                    // official: create new pose for each other pose in currentAssembly
                    // MVP impl: check if any given pose improves currentAssembly
                    bool swappedLigandPose = false;
                    for (const UniquePoseID &pose : ligand.getPoses()) {
                        // check if using this pose improves currentAssembly
                        if(pose.getLigandInternalPoseId() == currentAssembly.getPoseOfLigand(ligand.getID()))
                        {
                            continue;
                        }
                        LigandAlignmentAssembly assemblyCopy = currentAssembly;
                        assemblyCopy.swapPoseForLigand(ligand.getID(), pose.getLigandInternalPoseId());
                        //avoid identity swap
                        double const newAssemblyScore
                            = AssemblyScorer::calculateAssemblyScore(assemblyCopy, m_pairwiseAlignments, m_ligands);

                        if (newAssemblyScore > currentAssemblyScore) {
                            currentAssembly = assemblyCopy;
                            currentAssemblyScore = newAssemblyScore;
                            ligandAvailable.setAllAvailable();
                            swappedLigandPose = true;
                            break;
                        }
                    }
                    if(!swappedLigandPose)
                    {
                        ligandAvailable.at(ligand.getID()) = false;
                    }
                }

            }
            double const assemblyScore
                = AssemblyScorer::calculateAssemblyScore(currentAssembly, m_pairwiseAlignments, m_ligands);
            spdlog::info("Score after opt: {}", assemblyScore);
            if (assemblyScore > currentBestScore) {
                currentBestAssembly = currentAssembly;
                currentBestScore = assemblyScore;
            }
        }

        return {currentBestScore, currentBestAssembly.getAssemblyMapping(), m_ligands};
    }

}  // namespace coaler::multialign
