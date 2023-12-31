/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#include "Ligand.hpp"

namespace coaler::multialign {

    Ligand::Ligand(const RDKit::ROMol& mol, const UniquePoseSet& poses, LigandID id)
        : m_molecule(boost::make_shared<RDKit::ROMol>(mol)), m_poses(poses), m_id(id) {}

    /*----------------------------------------------------------------------------------------------------------------*/

    UniquePoseSet Ligand::getPoses() const noexcept { return m_poses; }

    /*----------------------------------------------------------------------------------------------------------------*/

    LigandID Ligand::getID() const noexcept { return m_id; }

    /*----------------------------------------------------------------------------------------------------------------*/

    unsigned Ligand::getNumHeavyAtoms() const noexcept { return m_molecule->getNumHeavyAtoms(); }

    /*----------------------------------------------------------------------------------------------------------------*/

    unsigned Ligand::getNumPoses() const noexcept { return m_poses.size(); }

    /*----------------------------------------------------------------------------------------------------------------*/

    RDKit::ROMol Ligand::getMolecule() const noexcept { return *m_molecule; }

    /*----------------------------------------------------------------------------------------------------------------*/

    RDKit::ROMOL_SPTR Ligand::getMoleculePtr() const noexcept { return m_molecule; }

    /*----------------------------------------------------------------------------------------------------------------*/

    void Ligand::addPose(const PoseID& poseId) noexcept { m_poses.insert(UniquePoseID(this->getID(),poseId)); }

}  // namespace coaler::multialign
