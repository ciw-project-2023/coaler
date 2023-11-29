//
// Created by chris on 11/5/23.
//

#include "MultiAligner.hpp"

#include <spdlog/spdlog.h>

#include <queue>
#include <utility>

#include "AssemblyScorer.hpp"
#include "LigandAlignmentAssembly.hpp"
#include "StartingAssemblyGenerator.hpp"

namespace {  // TODO move to own class?

    using LigandAvailabilityMapping = std::unordered_map<coaler::multialign::LigandID, bool>;

    void set_all_ligands_to_available(LigandAvailabilityMapping& mapping,
                                      const coaler::multialign::LigandVector& ligands)
    {
        mapping.clear();
        for(const coaler::multialign::Ligand& ligand : ligands)
        {
            mapping.emplace(ligand.getID(), true);
        }
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    /*auto LigandIsAvailableLambda = [](std::pair<coaler::multialign::LigandID, bool> entry){
      return entry.second;
    }; */

    /*----------------------------------------------------------------------------------------------------------------*/

    struct LigandIsAvailable{

        bool operator()(std::pair<coaler::multialign::LigandID, bool> entry)
        {
            return entry.second;
        }
    };

    /*----------------------------------------------------------------------------------------------------------------*/

    coaler::multialign::PairwiseAlignment calculate_alignment_scores(
        const std::vector<coaler::multialign::Ligand>& ligands,
        const std::shared_ptr<coaler::SingleAligner>& alignerPtr, const RDKit::ROMol& core) {
        coaler::multialign::PairwiseAlignment scores;
        for (coaler::multialign::LigandID firstMolId = 0; firstMolId < ligands.size();
             firstMolId++)  // TODO use ligand vector iterator
        {
            for (coaler::multialign::LigandID secondMolId = firstMolId + 1; secondMolId < ligands.size();
                 secondMolId++) {
                unsigned nofPosesFirst = ligands.at(firstMolId).getNofPoses();
                unsigned nofPosesSecond = ligands.at(secondMolId).getNofPoses();

                for (unsigned firstMolPoseId = 0; firstMolPoseId < nofPosesFirst; firstMolPoseId++) {
                    for (unsigned secondMolPoseId = 0; secondMolPoseId < nofPosesSecond; secondMolPoseId++) {
                        double alignmentScore = alignerPtr->align_molecules_kabsch(
                            ligands.at(firstMolId).getMolecule(), ligands.at(secondMolId).getMolecule(), firstMolPoseId,
                            secondMolPoseId, core);
                        coaler::multialign::UniquePoseIdentifier firstPose(firstMolId, firstMolPoseId);
                        coaler::multialign::UniquePoseIdentifier secondPose(secondMolId, secondMolPoseId);
                        scores.emplace(coaler::multialign::PosePair(firstPose, secondPose), alignmentScore);
                    }
                }
            }
        }
        return scores;
    }
}  // namespace

namespace coaler::multialign {

    using AssemblyWithScore = std::pair<LigandAlignmentAssembly, double>;

    struct AssemblyWithScoreLess {
        bool operator()(const AssemblyWithScore& lhs, const AssemblyWithScore& rhs) {
            if (lhs.first.getMissingLigandsCount() != lhs.first.getMissingLigandsCount()) {
                return lhs.first.getMissingLigandsCount() < lhs.first.getMissingLigandsCount();
            }
            return lhs.second < rhs.second;

        }
    };

    /*----------------------------------------------------------------------------------------------------------------*/

    struct AssemblyWithScoreGreater {
        bool operator()(const AssemblyWithScore& lhs, const AssemblyWithScore& rhs) {
            if (lhs.first.getMissingLigandsCount() != lhs.first.getMissingLigandsCount()) {
                return lhs.first.getMissingLigandsCount() > lhs.first.getMissingLigandsCount();
            }
            return lhs.second > rhs.second;

        }
    };

    /*----------------------------------------------------------------------------------------------------------------*/

    MultiAligner::MultiAligner(const std::vector<RDKit::RWMol>& molecules, RDKit::ROMol core,
                               const coaler::SingleAligner& aligner, unsigned maxStartingAssemblies)
        : m_core(std::move(core)), m_singleAligner(aligner), m_maxStartingAssemblies(maxStartingAssemblies) {
        // create ligands
        for (LigandID id = 0; id < molecules.size(); id++) {
            UniquePoseSet poses;
            for (PoseID poseId = 0; poseId < molecules.at(id).getNumConformers(); poseId++) {
                poses.emplace(id, poseId);
            }
            m_ligands.emplace_back(molecules.at(id), poses, id);
        }
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    MultiAlignerResult MultiAligner::alignMolecules() {
        // calculate pairwise alignments
        PairwiseAlignment allPosesAlignmentScores
            = calculate_alignment_scores(m_ligands, std::make_shared<coaler::SingleAligner>(m_singleAligner), m_core);

        m_poseRegisters = PoseRegisterBuilder::buildPoseRegisters(allPosesAlignmentScores, m_ligands);

        // build starting ensembles from registers
        // AssemblyCollection assemblies;
        std::priority_queue<AssemblyWithScore, std::vector<AssemblyWithScore>, AssemblyWithScoreLess> assemblies;

        for (const Ligand& ligand : m_ligands) {
            for (const UniquePoseIdentifier& pose : ligand.getPoses()) {
                LigandAlignmentAssembly assembly
                    = StartingAssemblyGenerator::generateStartingAssembly(pose, m_poseRegisters, m_ligands);
                AssemblyWithScore assemblyWithScore = std::make_pair(
                    assembly, AssemblyScorer::calculateAssemblyScore(assembly, m_pairwiseAlignments, m_ligands));
                // insert if queue no full or new assembly is larger that worst assembly in queue
                if (assemblies.size() < m_maxStartingAssemblies || AssemblyWithScoreLess()(assemblies.top(), assemblyWithScore)) { //TODO ensure this is called correctly
                    assemblies.push(assemblyWithScore);
                }
            }
        }
        // top #m_maxStartingAssemblies are now in queue. find best scoring assembly by optimizing all

        //optimize all starting assemblies.

        LigandAlignmentAssembly currentBestAssembly = assemblies.top().first; //TODO default constructor for assembly?
        double currentBestScore = AssemblyScorer::calculateAssemblyScore(currentBestAssembly,
                                                                  m_pairwiseAlignments,
                                                                  m_ligands);
        while(!assemblies.empty())
        {
            LigandAlignmentAssembly assembly = assemblies.top().first;
            assemblies.pop();
            LigandAvailabilityMapping ligandAvailabilities; // TODO init here and only set all true in anon function
            set_all_ligands_to_available(ligandAvailabilities, m_ligands);
            while(std::any_of(ligandAvailabilities.begin(), ligandAvailabilities.end(), LigandIsAvailable()))
            {
                for(const Ligand& ligand : m_ligands)
                {
                    if(!ligandAvailabilities.at(ligand.getID()))
                    {
                        continue;
                    }
                    // official: create new pose for each other pose in assembly

                    //MVP impl: check if any given pose improves assembly
                    for(const UniquePoseIdentifier& pose : ligand.getPoses())
                    {
                        //check if using this pose improves assembly
                        LigandAlignmentAssembly assemblyCopy = assembly;
                        assemblyCopy.swapPoseForLigand(ligand.getID(), pose.getLigandInternalPoseId());
                        double newAssemblyScore = AssemblyScorer::calculateAssemblyScore(assemblyCopy,
                                                                                         m_pairwiseAlignments,
                                                                                         m_ligands);
                        double currentAssemblyScore = AssemblyScorer::calculateAssemblyScore(assembly,
                                                                                            m_pairwiseAlignments,
                                                                                            m_ligands);
                        if(newAssemblyScore > currentAssemblyScore)
                        {
                            assembly = assemblyCopy;
                            set_all_ligands_to_available(ligandAvailabilities, m_ligands);
                        }
                    }

                }
            }
            double assemblyScore = AssemblyScorer::calculateAssemblyScore(assembly,
                                                                         m_pairwiseAlignments,
                                                                         m_ligands);
            if(assemblyScore > currentBestScore)
            {
                currentBestAssembly = assembly;
                currentBestScore = assemblyScore;
            }
        }

        return MultiAlignerResult(currentBestScore, currentBestAssembly.getAssemblyMapping(), m_ligands);
    }

}  // namespace coaler::multialign