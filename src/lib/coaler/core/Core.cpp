//
// Created by malte on 12/4/23.
//

#include "Core.hpp"
#include "GraphMol/FMCS/FMCS.h"
#include "GraphMol/Substruct/SubstructMatch.h"
#include "GraphMol/ChemTransforms/ChemTransforms.h"

namespace coaler::core {


    Core::Core(const std::vector<RDKit::RWMol>& molecules, const coreType coreType)
        :m_molecules(std::move(molecules)), m_coreType(coreType) {
        assert(m_molecules.size() > 0);
    }

    coreType Core::getScaffoldType() const {
        return m_coreType;
    }

    CoreAsMol Core::getCore() const {
        return m_core;
    }

    std::vector<RDKit::RWMol> Core::calculateCoreMcs() const {

        std::vector<RDKit::ROMOL_SPTR> shared_mols;
        for(RDKit::RWMol mol : m_molecules) {
            shared_mols.emplace_back(boost::make_shared<RDKit::ROMol>(mol));
        }


        RDKit::MCSResult mcs = RDKit::findMCS(shared_mols);
        RDKit::ROMOL_SPTR mcs_structure = mcs.QueryMol;

        std::vector<RDKit::RWMol> result;

        for(RDKit::RWMol mol : m_molecules) {
            auto match = RDKit::SubstructMatch(mol, *mcs_structure);
        }
    }

    std::vector<RDKit::RWMol> Core::calculateCoreMurcko() const {

        std::vector<RDKit::ROMOL_SPTR> shared_mols;
        for(RDKit::RWMol mol : m_molecules) {
            shared_mols.emplace_back(boost::make_shared<RDKit::ROMol>(mol));
        }

        RDKit::MCSResult mcs = RDKit::findMCS(shared_mols);
        RDKit::ROMOL_SPTR mcs_structure = mcs.QueryMol;


        std::vector<RDKit::RWMol> result;

        for(RDKit::RWMol mol : m_molecules) {
            RDKit::ROMol* scaffold = RDKit::MurckoDecompose(mol);

            result.emplace_back(*scaffold);
        }

        return result;
    }
