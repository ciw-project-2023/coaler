#include <GraphMol/SmilesParse/SmartsWrite.h>

#include <cstdint>
#include <filesystem>

#include "../src/parser/FileParser.hpp"
#include "../src/singlealign/SingleAligner.hpp"
#include "catch2/catch.hpp"

TEST_CASE("Single Aligner", "[aligner]") {
    SECTION("Error when core structure is too small!") {
        RDKit::RWMol *mol_a = RDKit::SmilesToMol("CCCN");
        RDKit::RWMol *mol_b = RDKit::SmilesToMol("CCCO");

        coaler::SingleAligner single_aligner(5, 80);
        CHECK_THROWS_WITH(single_aligner.align_molecules_kabsch(*mol_a, *mol_b, std::nullopt), "Size of core is too small!");
    };
    SECTION("Warning when core structure is too large!") {
        RDKit::RWMol *mol_a = RDKit::SmilesToMol("CCCN");
        RDKit::RWMol *mol_b = RDKit::SmilesToMol("CCCO");

        coaler::SingleAligner single_aligner(1, 1);
        single_aligner.align_molecules_kabsch(*mol_a, *mol_b, std::nullopt);
    };
    SECTION("Core is only part of one molecule") {
        RDKit::RWMol *mol_a = RDKit::SmilesToMol("CCCN");
        RDKit::RWMol *mol_b = RDKit::SmilesToMol("CCCO");

        RDKit::RWMol *core = RDKit::SmilesToMol("CO");

        coaler::SingleAligner single_aligner(1, 80);
        CHECK_THROWS_WITH(single_aligner.align_molecules_kabsch(*mol_a, *mol_b, *core),
                          "Core is not a common core structure of molecule a and molecule b!");
    };
    SECTION("Compare two molecules") {
        RDKit::RWMol *mol_a = RDKit::SmilesToMol("CCCN");
        RDKit::RWMol *mol_b = RDKit::SmilesToMol("CCCO");

        coaler::SingleAligner single_aligner(1, 80);
        auto result = single_aligner.align_molecules_kabsch(*mol_a, *mol_b, std::nullopt);

        double score = static_cast<unsigned int>((std::get<0>(result) * 10000)) / 10000.0;
        CHECK(score == 0.0127);

        RDKit::ROMOL_SPTR core = std::get<1>(result);
        CHECK(RDKit::MolToSmarts(*core) == "[#6]-[#6]-[#6]");
    };
}
