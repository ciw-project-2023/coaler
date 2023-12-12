/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#include "Ligand.hpp"

namespace coaler::multialign {

    Ligand::Ligand(const RDKit::RWMol& mol, const UniquePoseSet& poses, LigandID id)
        : m_molecule(mol), m_poses(poses), m_id(id) {}

    /*----------------------------------------------------------------------------------------------------------------*/

    UniquePoseSet Ligand::getPoses() const noexcept { return m_poses; }

    /*----------------------------------------------------------------------------------------------------------------*/

    LigandID Ligand::getID() const noexcept { return m_id; }

    /*----------------------------------------------------------------------------------------------------------------*/

    unsigned Ligand::getNumHeavyAtoms() const noexcept { return m_molecule.getNumHeavyAtoms(); }

    /*----------------------------------------------------------------------------------------------------------------*/

    unsigned Ligand::getNumPoses() const noexcept { return m_poses.size(); }

    RDKit::RWMol Ligand::getMolecule() const noexcept { return m_molecule; }

    void Ligand::setMolecule(RDKit::RWMol new_mol) {
        m_molecule = new_mol;
    }

}  // namespace coaler::multialign
