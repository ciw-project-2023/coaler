#include "SingleAligner.hpp"

#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/FMCS/FMCS.h>
#include <RDGeneral/export.h>
#include <vector>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>

namespace ciw {

    SingleAligner::SingleAligner() {
    }

    SingleAligner::~SingleAligner() {
    }

    void SingleAligner::set_outputfile(std::string) {

    }

    void SingleAligner::algin_molecules(RDKit::ROMol mol_a, RDKit::ROMol mol_b) {
        // zip mol_a and mol_b into a vector
        std::vector<RDKit::ROMOL_SPTR> mols;
        mols.push_back(RDKit::ROMOL_SPTR(new RDKit::ROMol(mol_a)));
        mols.push_back(RDKit::ROMOL_SPTR(new RDKit::ROMol(mol_b)));

        // find the MCS
        RDKit::MCSResult res = RDKit::findMCS(mols);

        // print the MCS
        std::cout << "MCS: " << res.SmartsString << std::endl;
    }

} // ciw