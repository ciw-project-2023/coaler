//
// Created by chris on 11/5/23.
//

#include "MultiAligner.hpp"

#include <utility>
#include <spdlog/spdlog.h>

#include "LigandAlignmentAssembly.hpp"
#include "StartingAssemblyGenerator.hpp"

namespace{ //TODO move to own class?
    MultiAlign::PairwiseAlignment calculate_alignment_scores(
            const std::vector<MultiAlign::Ligand>& ligands,
            const std::shared_ptr<coaler::SingleAligner>& alignerPtr,
            const RDKit::ROMol& core)
    {
        MultiAlign::PairwiseAlignment scores;
        for(MultiAlign::LigandID firstMolId = 0; firstMolId < ligands.size(); firstMolId++) //TODO use ligand vector iterator
        {
            for(MultiAlign::LigandID secondMolId = firstMolId + 1; secondMolId < ligands.size(); secondMolId++)
            {
                unsigned nofPosesFirst = ligands.at(firstMolId).getNofPoses();
                unsigned nofPosesSecond = ligands.at(secondMolId).getNofPoses();

                 for(unsigned firstMolPoseId = 0; firstMolPoseId < nofPosesFirst; firstMolPoseId++)
                 {
                     for(unsigned  secondMolPoseId = 0; secondMolPoseId < nofPosesSecond; secondMolPoseId++)
                     {
                         double alignmentScore = alignerPtr->align_molecules_kabsch(
                                 ligands.at(firstMolId).getMolecule(),
                                 ligands.at(secondMolId).getMolecule(),
                                 firstMolPoseId,
                                 secondMolPoseId,
                                 core);
                         MultiAlign::UniquePoseIdentifier firstPose(firstMolId, firstMolPoseId);
                         MultiAlign::UniquePoseIdentifier secondPose(secondMolId, secondMolPoseId);
                         scores.emplace(MultiAlign::PosePair(firstPose, secondPose), alignmentScore);
                     }
                 }
            }
        }
        return scores;
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    double calc_score_of_assembly( //TODO move to assemlby and store to avoid recalcing a lot?
        const MultiAlign::LigandAlignmentAssembly& assembly,
        const MultiAlign::PairwiseAlignment& scores,
        const MultiAlign::LigandVector& ligands)
    {
        double assemblyScore = 0.0;
        for(const MultiAlign::Ligand& firstLigand : ligands)
        {
            for(const MultiAlign::Ligand& secondLigand : ligands)
            {
                if(firstLigand.getID() >= secondLigand.getID())
                {
                    continue;
                }
                MultiAlign::PoseID firstLigandPoseID = assembly.getPoseOfLigand(firstLigand.getID());
                MultiAlign::PoseID secondLigandPoseID = assembly.getPoseOfLigand(secondLigand.getID());

                if(firstLigandPoseID == std::numeric_limits<MultiAlign::PoseID>::max()
                    || secondLigandPoseID == std::numeric_limits<MultiAlign::PoseID>::max())
                {
                    spdlog::info("encountered invalid PosePair during optimization.");
                    continue;
                }

                MultiAlign::UniquePoseIdentifier firstLigandPose{
                    firstLigand.getID(), firstLigandPoseID};
                MultiAlign::UniquePoseIdentifier secondLigandPose{
                    secondLigand.getID(), secondLigandPoseID};
                assemblyScore += scores.at({firstLigandPose, secondLigandPose});
            }
        }
        return assemblyScore;
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    double get_score_deficit_of_ligand(
        const MultiAlign::Ligand& ligand,
        const MultiAlign::LigandVector& ligands,
        const MultiAlign::LigandAlignmentAssembly& assembly,
        const MultiAlign::PairwiseAlignment& scores)
    {
        double scoreDeficit = 0.0;
        for(const MultiAlign::Ligand& otherLigand : ligands)
        {
            if(ligand.getID() == otherLigand.getID())
            {
                continue;
            }
            MultiAlign::UniquePoseIdentifier ligandPose = {ligand.getID(),
            assembly.getPoseOfLigand(ligand.getID())};
        }
    }
}

namespace MultiAlign
{

    using AssemblyCollection = std::unordered_map<UniquePoseIdentifier, LigandAlignmentAssembly, UniquePoseIdentifierHash>;

    MultiAligner::MultiAligner(
            const std::vector<RDKit::RWMol>& molecules,
            RDKit::ROMol  core,
            const coaler::SingleAligner& aligner)
            : m_core(std::move(core))
            , m_singleAligner(aligner)
            {
                //create ligands
                for(LigandID id = 0; id < molecules.size(); id++)
                {
                    UniquePoseSet poses;
                    for(PoseID poseId = 0; poseId < molecules.at(id).getNumConformers(); poseId++)
                    {
                        poses.emplace(id, poseId);
                    }
                    m_ligands.emplace_back(molecules.at(id),
                                         poses,
                                         id);
                }
            }

    /*----------------------------------------------------------------------------------------------------------------*/

    MultiAlignerResult MultiAligner::alignMolecules(){
        //calculate pairwise alignments
        PairwiseAlignment allPosesAlignmentScores = calculate_alignment_scores(
                m_ligands,
                std::make_shared<coaler::SingleAligner>(m_singleAligner),
                m_core);

        m_poseRegisters = PoseRegisterBuilder::buildPoseRegisters(
                allPosesAlignmentScores,
                m_ligands);

        // build starting ensembles from registers
        AssemblyCollection assemblies;
        for(const Ligand& ligand : m_ligands)
        {
            for(const UniquePoseIdentifier& pose : ligand.getPoses())
            {
                assemblies.emplace(pose, StartingAssemblyGenerator::generateStartingAssembly(
                    pose,
                    m_poseRegisters,
                    m_ligands
                    ));
            }
        }

        //optimize ensembles

        /**
         * forall assemblies
         *      while [not aboirt condition]
         *          get ligand with worst score
         *                 forall other ligands
         *                 generate new pose of ligand aligned to poses of other ligand
         *                 find new ligand pose with best fit
         *                 exchange pose in assembly
         *
         *
         *
         *
         */

        return MultiAlignerResult();
    }
}