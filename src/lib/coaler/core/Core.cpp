//
// Created by malte on 12/4/23.
//

#include "Core.hpp"
#include "GraphMol/FMCS/FMCS.h"
#include "GraphMol/ChemTransforms/ChemTransforms.h"
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>

namespace coaler::core {

    Core::Core(const std::vector<RDKit::RWMol>& molecules, const coreType coreType)
        : m_molecules(molecules), m_coreType(coreType) {
        assert(!m_molecules.empty());

        if (coreType == coreType::MCS) {
            m_core = calculateCoreMcs();
        }
        else {
            m_core = calculateCoreMurcko();
        }

    }


    coreType Core::getScaffoldType() const { return m_coreType; }


    CoreAsMol Core::getCore() const { return m_core; }

    CoreAsMol Core::calculateCoreMcs() const {
        std::vector<RDKit::ROMOL_SPTR> shared_mols;
        for (RDKit::RWMol mol : m_molecules) {
            shared_mols.emplace_back(boost::make_shared<RDKit::ROMol>(mol));
        }

        RDKit::MCSResult const mcs = RDKit::findMCS(shared_mols);
        RDKit::ROMOL_SPTR const mcs_structure = mcs.QueryMol;

        // Ja das hier ist extrem hässlich aber ich habe noch keine bessere Möglichkeit gefunden
        // Die Methode findMCS() scheint bei der Erstellung vom ROMOL keine Infos über die Ringstruktur
        // zu machen und daher gibt MurkoDecompose Errors aus.
        const std::string molAsSmarts = RDKit::MolToSmarts(*mcs_structure);
        RDKit::ROMol molAsROMol = *RDKit::SmilesToMol(molAsSmarts);

        return molAsROMol;
    }

    CoreAsMol Core::calculateCoreMurcko() const {
        RDKit::ROMol const mcs_structure = calculateCoreMcs();

        //CoreAsMol ret = *RDKit::MurckoDecompose(*RDKit::SmilesToMol("c1ccncc1"));
        CoreAsMol ret = *RDKit::MurckoDecompose(mcs_structure);
        return ret;
    }


}