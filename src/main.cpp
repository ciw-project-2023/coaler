#include <cstdint> // todo: is necessary when using RDKit libraries
#include <iostream>

#include <spdlog/spdlog.h>

#include <GraphMol/SmilesParse/SmilesParse.h>

#include "SingleAligner.hpp"

int main(int argc, char *argv[]) {
    RDKit::RWMol *mol_a = RDKit::SmilesToMol("CCCO");
    RDKit::RWMol *mol_b = RDKit::SmilesToMol("CCCN");

    spdlog::info("Starting program");

    spdlog::info("No. of Atoms: {}", mol_a->getNumAtoms());

    ciw::SingleAligner single_aligner;

    spdlog::info("Program linked and executed correctly");

}
