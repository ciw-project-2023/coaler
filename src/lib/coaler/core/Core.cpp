//
// Created by malte on 12/4/23.
//

#include "Core.hpp"

#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include "GraphMol/ChemTransforms/ChemTransforms.h"
#include "GraphMol/FMCS/FMCS.h"

namespace coaler::core {

    Core::Core() {}

    Core::Core(const std::vector<RDKit::RWMol>& molecules, const CoreType coreType)
        : m_molecules(molecules), m_coreType(coreType) {
        assert(!m_molecules.empty());

        if (coreType == CoreType::MCS) {
            m_core = calculateCoreMcs();
        } else {
            m_core = calculateCoreMurcko();
        }
    }

    CoreType Core::getCoreType() const { return m_coreType; }

    CoreAsMol Core::getCore() const { return m_core; }

    CoreAsMol Core::calculateCoreMcs() const {
        std::vector<RDKit::ROMOL_SPTR> shared_mols;
        for (const RDKit::RWMol& mol : m_molecules) {
            shared_mols.emplace_back(boost::make_shared<RDKit::ROMol>(mol));
        }

        RDKit::MCSResult const mcs = RDKit::findMCS(shared_mols);
        RDKit::RWMol mcsAsMol = *mcs.QueryMol;

        //RDKit::MolOps::sanitizeMol(mcsAsMol);

        return mcsAsMol;
    }

    CoreAsMol Core::calculateCoreMurcko() const {
        RDKit::ROMol const mcs_structure = calculateCoreMcs();

        // CoreAsMol ret = *RDKit::MurckoDecompose(*RDKit::SmilesToMol("c1ccncc1"));
        CoreAsMol ret = *RDKit::MurckoDecompose(mcs_structure);
        ret.RDKit::ROMol::clearComputedProps();
        ret.RDKit::ROMol::updatePropertyCache();
        return ret;
    }

}  // namespace coaler::core