//
// Created by chris on 11/4/23.
//
#pragma once

#include "Forward.hpp"
#include "UniquePoseIdentifier.hpp"
#include <GraphMol/GraphMol.h>

#include <unordered_set>

namespace MultiAlign {

    class Ligand {
    public:
        Ligand(const RDKit::RWMol& mol,
               const UniquePoseSet &poses,
               LigandID id);

        [[nodiscard]] UniquePoseSet getPoses() const noexcept;

        [[nodiscard]] LigandID getID() const noexcept;

        unsigned getHeavyAtomSize() const noexcept;

    private:
        LigandID m_id;
        RDKit::RWMol m_molecule;
        UniquePoseSet m_poses;
    };

}
