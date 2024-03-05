#pragma once

#include <GraphMol/GraphMol.h>

#include "Alias.hpp"
#include "UniquePoseID.hpp"
#include "UniquePoseSet.hpp"
/*!
 * @file Ligand.hpp
 * @brief This file contains the Ligand class which is used to represent a ligand molecule and its conformers.
 */
namespace coaler::multialign {
    /**
     * The Ligand class provides functionality for the representation of a ligand molecule and its conformers.
     */
    class Ligand {
      public:
        Ligand(const RDKit::ROMol& mol, const UniquePoseSet& poses, LigandID id);

        /**
         * get idenitifers of all poses embedded in ligand.
         * @return The ids of the ligands molecule conformers.
         */
        [[nodiscard]] UniquePoseSet getPoses() const noexcept;

        /**
         * Add a new pose to the map
         * @param poseId The ids of the ligands molecule conformer.
         */
        void addPose(const PoseID& poseId) noexcept;

        /**
         * Get the id of the ligand.
         * @return The id of the ligand.
         */
        [[nodiscard]] LigandID getID() const noexcept;

        /**
         * Get the number of heavy atoms in the ligand.
         * @return The number of heavy atoms in the ligand.
         */
        [[nodiscard]] unsigned getNumHeavyAtoms() const noexcept;

        /**
         * Get the number of conformers embedded in the ligand.
         * @return The number of conformers embedded in the ligand.
         */
        [[nodiscard]] unsigned getNumPoses() const noexcept;

        /**
         *
         * @return The molecule represented by the ligand.
         */
        [[nodiscard]] RDKit::ROMol getMolecule() const noexcept;

        [[nodiscard]] RDKit::RWMol const* getMoleculePtr() const noexcept;

        void removePose(PoseID pose);

        bool operator==(const Ligand& other) const;

      private:
        LigandID m_id;
        RDKit::RWMol m_molecule;
        UniquePoseSet m_poses;
    };

}  // namespace coaler::multialign
