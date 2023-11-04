#include <cstdint> // todo: is necessary when using RDKit libraries
#include <iostream>

#include <GraphMol/SmilesParse/SmilesParse.h>

#include "SingleAligner.hpp"

int main(int argc, char *argv[]) {
    //Using two cox2 inhibitor molecules as example
    //Celebrex
    RDKit::RWMol *mol_a = RDKit::SmilesToMol("CC1=CC=C(C=C1)C1=CC(=NN1C1=CC=C(C=C1)S(N)(=O)=O)C(F)(F)F");
    //Bextra
    RDKit::RWMol *mol_b = RDKit::SmilesToMol("CC1=C(C(=NO1)C1=CC=CC=C1)C1=CC=C(C=C1)S(N)(=O)=O");

    // print number of atoms
    std::cout << "Number of atoms in Celebrex: " << mol_a->getNumAtoms() << std::endl;
    std::cout << "Number of atoms in Bextra: " << mol_b->getNumAtoms() << std::endl;
    ciw::SingleAligner single_aligner;

    // align molecules
    single_aligner.align_molecules_kabsch(*mol_a, *mol_b, std::nullopt);
}