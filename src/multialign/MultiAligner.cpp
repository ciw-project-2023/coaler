//
// Created by chris on 11/5/23.
//

#include "MultiAligner.hpp"

namespace{
    MultiAlign::PairwiseAlignment calculate_alignment_scores(
            const std::vector<RDKit::RWMol>& molecules,
            std::shared_ptr<coaler::SingleAligner> alignerPtr,
            const RDKit::ROMol& core)
    {
        MultiAlign::PairwiseAlignment scores;
        for(MultiAlign::LigandID firstMolId = 0; firstMolId < molecules.size(); firstMolId++)
        {
            for(MultiAlign::LigandID secondMolId = firstMolId + 1; secondMolId < molecules.size(); secondMolId++)
            {
                unsigned nofPosesFirst = molecules.at(firstMolId).getNumConformers();
                unsigned nofPosesSecond = molecules.at(secondMolId).getNumConformers();

                 for(unsigned firstMolPoseId = 0; firstMolPoseId < nofPosesFirst; firstMolPoseId++)
                 {
                     for(unsigned  secondMolPoseId = 0; secondMolPoseId < nofPosesSecond; secondMolPoseId++)
                     {
                         double alignmentScore = alignerPtr->align_molecules_kabsch(
                                 molecules.at(firstMolId),
                                 molecules.at(secondMolId),
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
    }
}

namespace MultiAlign
{
    MultiAligner::MultiAligner(
            const std::vector<RDKit::RWMol>& molecules,
            const RDKit::ROMol& core,
            const coaler::SingleAligner& aligner)
            : m_molecules(molecules)
            , m_core(core)
            , m_singleAligner(aligner)
            {
            }

    /*----------------------------------------------------------------------------------------------------------------*/

    MultiAlignerResult MultiAligner::alignMolecules(){
        PairwiseAlignment allPosesAlignmentScores = calculate_alignment_scores(
                m_molecules,
                std::make_shared<coaler::SingleAligner>(m_singleAligner),
                m_core);

        std::vector<Ligand> ligands;
        for(LigandID id = 0; id < m_molecules.size(); id++)
        {
            UniquePoseSet poses;
            for(PoseID poseId = 0; poseId < m_molecules.at(id).getNumConformers(); poseId++)
            {
                poses.emplace(id, poseId);
            }
            ligands.emplace_back(m_molecules.at(id),
                            poses,
                            id);
        }

        m_poseRegisters = PoseRegisterBuilder::buildPoseRegisters(
                allPosesAlignmentScores,
                ligands);
        assert(!m_poseRegisters.empty());


        //generate registers

        // build starting ensembles from registers

        //optimize ensembles

        return MultiAlignerResult();
    }
}