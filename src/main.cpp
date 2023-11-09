#include <cstdint> // todo: is necessary when using RDKit libraries

#include <spdlog/spdlog.h>

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/MolSupplier.h>

#include "SingleAligner.hpp"

int main(int argc, char *argv[]) {
    RDKit::RWMol *mol_a = RDKit::SmilesToMol("CCCO");
    RDKit::RWMol *mol_b = RDKit::SmilesToMol("CCCN");

    // Read the molecule from the SDF file.
    //RDKit::SDMolSupplier supplier(input_sdf_file_path);
    //RDKit::ROMol* molecule = supplier.next();

    spdlog::info("Starting program");

    spdlog::info("No. of Atoms: {}", mol_a->getNumAtoms());

    ciw::SingleAligner single_aligner(5, 50);

    spdlog::info("Program linked and executed correctly");

}
