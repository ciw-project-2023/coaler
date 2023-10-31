#include "SingleAligner.hpp"

#include <GraphMol/MolAlign/AlignMolecules.h>

namespace ciw {

    SingleAligner::SingleAligner() {
    }

    SingleAligner::~SingleAligner() {
    }

    void SingleAligner::set_outputfile(std::string) {

    }

    void SingleAligner::algin_molecules(RDKit::ROMol mol_a, RDKit::ROMol mol_b) {
        //RDKit::MolAlign::getAlignmentTransform(mol_a, mol_b);
    }

} // ciw