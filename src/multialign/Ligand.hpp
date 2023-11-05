//
// Created by chris on 11/4/23.
//
#pragma once

#include "Forward.hpp"

#include <set>

namespace MultiAlign {

    class Ligand {
    public:
        Ligand(const std::set<PoseID> &poses,
               LigandID id);

        [[nodiscard]] std::set<PoseID> getPoses() const noexcept;

        [[nodiscard]] LigandID getID() const noexcept;

    private:
        LigandID m_id;

        std::set<PoseID> m_poses;

    };

}
