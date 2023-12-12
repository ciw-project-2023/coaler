/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */
#pragma once

#include <GraphMol/GraphMol.h>

#include <unordered_set>

#include "Forward.hpp"
#include "UniquePoseID.hpp"

namespace coaler::multialign {

    /**
     * Represents a molecule.
     */
    class Ligand {
      public:
        Ligand(const RDKit::RWMol& mol, const UniquePoseSet& poses, LigandID id);

        /**
         * get idenitifers of all poses embedded in ligand.
         * @return
         */
        [[nodiscard]] UniquePoseSet getPoses() const noexcept;

        /**
         *
         * @return
         */
        [[nodiscard]] LigandID getID() const noexcept;

        /**
         *
         * @return
         */
        unsigned getNumHeavyAtoms() const noexcept;

        /**
         *
         * @return The number of conformers embedded in the ligand.
         */
        unsigned getNumPoses() const noexcept;

        /**
         *
         * @return The molecule represented by the ligand.
         */
        RDKit::RWMol getMolecule() const noexcept;

        void setMolecule(RDKit::RWMol new_mol);

      private:
        LigandID m_id;
        RDKit::RWMol m_molecule;
        UniquePoseSet m_poses;
    };

}  // namespace coaler::multialign
