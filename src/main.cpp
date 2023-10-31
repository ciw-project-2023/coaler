#include <cstdint> // todo: is necessary when using RDKit libraries
#include <iostream>

#include <GraphMol/SmilesParse/SmilesParse.h>

#include "SingleAligner.hpp"

int main(int argc, char *argv[]) {
    RDKit::RWMol* mol_a = RDKit::SmilesToMol("CCCO");
    RDKit::RWMol* mol_b = RDKit::SmilesToMol("CCCN");

    std::cout << mol_a->getNumAtoms() << std::endl;

    ciw::SingleAligner single_aligner;

    //single_aligner.algin_molecules(mol_a->, mol_b);

    std::cout << "RDKit is linked :)" << std::endl;
}