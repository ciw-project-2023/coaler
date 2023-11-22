#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <coaler/io/OutputWriter.hpp>

#include "catch2/catch.hpp"

TEST_CASE("Output Parser", "[io]") {
    coaler::io::OutputWriter output_parser;

    SECTION("Add one pair to output") {
        RDKit::RWMol *mol_a = RDKit::SmilesToMol("CCCN");
        RDKit::RWMol *mol_b = RDKit::SmilesToMol("CCCO");

        RDKit::DGeomHelpers::EmbedMolecule(*mol_a);
        RDKit::DGeomHelpers::EmbedMolecule(*mol_b);

        coaler::io::AlignedMolPair mol_pair;
        mol_pair.mol_a = RDKit::ROMol{*mol_a};
        mol_pair.mol_b = RDKit::ROMol{*mol_b};
        mol_pair.id_mol_a = 0;
        mol_pair.id_mol_b = 1;
        mol_pair.align_score = 0.55;

        output_parser.add_aligned_mols(mol_pair);

        SECTION("save output in file") { output_parser.save_molecules_w_scores_in_file("output"); };
        SECTION("print output in log") { output_parser.print_in_log_molecules_w_scores(); };
    };
}
