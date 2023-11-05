//
// Created by chris on 11/4/23.
//

#include "Ligand.hpp"

namespace MultiAlign {

    Ligand::Ligand(const std::set<PoseID> &poses, LigandID id)
            : m_poses(poses), m_id(id) {
    }

    std::set<PoseID> Ligand::getPoses() const noexcept {
        return m_poses;
    }

    LigandID Ligand::getID() const noexcept {
        return m_id;
    }

}