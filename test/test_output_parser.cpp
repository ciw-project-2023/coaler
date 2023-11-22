#include <cstdint>

#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include "../src/output/OutputParser.hpp"
#include "catch2/catch.hpp"

TEST_CASE("Output Parser", "[output]") {
    coaler::OutputParser output_parser;

    SECTION("Add one pair to output") {
        RDKit::RWMol *mol_a = RDKit::SmilesToMol("CCCN");
        RDKit::RWMol *mol_b = RDKit::SmilesToMol("CCCO");

        RDKit::DGeomHelpers::EmbedMolecule(*mol_a);
        RDKit::DGeomHelpers::EmbedMolecule(*mol_b);

        coaler::AlignedMolPair mol_pair;
        mol_pair.mol_a = RDKit::ROMol{*mol_a};
        mol_pair.mol_b = RDKit::ROMol{*mol_b};
        mol_pair.align_score = 0.55;

        output_parser.add_aligned_mols(mol_pair);
    };

    output_parser.save_molecules_w_scores_in_file("output");
}
