#include <cstdint> // todo: is necessary when using RDKit libraries

#include <iostream>

#include "SingleAligner.hpp"

int main(int argc, char *argv[]) {

    RDKit::ROMol mol_a;
    RDKit::ROMol mol_b;

    ciw::SingleAligner single_aligner;

    single_aligner.algin_molecules(mol_a, mol_b);

    std::cout << "RDKit is linked :)" << std::endl;
}