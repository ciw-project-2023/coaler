/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#include "Ligand.hpp"

namespace coaler::multialign {

    Ligand::Ligand(const RDKit::ROMol& mol, const UniquePoseSet& poses, LigandID id)
        : m_molecule(mol), m_poses(poses), m_id(id) {}

    /*----------------------------------------------------------------------------------------------------------------*/

    UniquePoseSet Ligand::getPoses() const noexcept { return m_poses; }

    /*----------------------------------------------------------------------------------------------------------------*/

    LigandID Ligand::getID() const noexcept { return m_id; }

    /*----------------------------------------------------------------------------------------------------------------*/

    unsigned Ligand::getNumHeavyAtoms() const noexcept { return m_molecule.getNumHeavyAtoms(); }

    /*----------------------------------------------------------------------------------------------------------------*/

    unsigned Ligand::getNumPoses() const noexcept { return m_poses.size(); }

    /*----------------------------------------------------------------------------------------------------------------*/

    RDKit::ROMol Ligand::getMolecule() const noexcept { return m_molecule; }

    /*----------------------------------------------------------------------------------------------------------------*/

    RDKit::RWMol const* Ligand::getMoleculePtr() const noexcept { return &m_molecule; }

    /*----------------------------------------------------------------------------------------------------------------*/

    void Ligand::addPose(const PoseID& poseId) noexcept { m_poses.insert(UniquePoseID(this->getID(), poseId)); }

    /*----------------------------------------------------------------------------------------------------------------*/

    void Ligand::removePose(const PoseID pose) {
        m_molecule.removeConformer(pose);
        m_poses.erase({this->getID(), pose});
    }
    bool Ligand::operator==(const Ligand& other) const {
        return this->getID() == other.getID();
    }

}  // namespace coaler::multialign
