//
// Created by chris on 11/5/23.
//

#include "LigandAlignmentAssembly.hpp"


namespace MultiAlign {


    LigandAlignmentAssembly::LigandAlignmentAssembly(
            const std::unordered_map<LigandID, PoseID>& initialAssembly)
        : m_assembly(initialAssembly)
    {

    }

    /*----------------------------------------------------------------------------------------------------------------*/

    void LigandAlignmentAssembly::swapPoseForLigand(
            const LigandID ligandId,
            const PoseID newPoseId) {
        assert(m_assembly.count(ligandId) != 0);
        m_assembly.at(ligandId) = newPoseId;
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    PoseID LigandAlignmentAssembly::getPoseOfLigand(LigandID ligandId) {
        assert(m_assembly.count(ligandId) != 0);
        return m_assembly.at(ligandId);
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    void LigandAlignmentAssembly::incrementMissingLigandsCount() {
        m_missingLigandsCount++;
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    bool LigandAlignmentAssembly::insertLigandPose(LigandID ligand, PoseID pose) {
        return m_assembly.emplace(ligand, pose).second;
    }

}