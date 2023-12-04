//
// Created by malte on 12/4/23.
//

#include "Core.hpp"
#include "GraphMol/FMCS/FMCS.h"
// TODO: Add Chem lib Chem/Scaffolds/MurckoScaffold

namespace coaler::core {


    Core::Core(const std::vector<RDKit::RWMol>& molecules, const coreType coreType)
        :m_molecules(std::move(molecules)), m_coreType(coreType) {
        assert(m_molecules.size() > 0);

        bool isCalculated = calculateCore();
        assert(isCalculated);
    }

    coreType Core::getScaffoldType() const {
        return m_coreType;
    }

    CoreAsMol Core::getCore() const {
        return m_core;
    }

    bool Core::calculateCore() const {

        std::vector<RDKit::ROMOL_SPTR> shared_mols;
        for(RDKit::RWMol mol : m_molecules) {
            shared_mols.emplace_back(boost::make_shared<RDKit::ROMol>(mol));
        }


        RDKit::MCSResult mcs = RDKit::findMCS(shared_mols);
        RDKit::ROMOL_SPTR mcs_structure = mcs.QueryMol;

        if (m_coreType == coreType::MCS) {
            return true;
        }
        else if (m_coreType == coreType::Murcko) {
            RDKit::RWMol murcko = RDKit::getScaffoldForMol(mcs);
        }
    }

}