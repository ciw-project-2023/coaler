#include "catch2/catch.hpp"

#include "../src/SingleAligner.hpp"
#include "../src/parsers/FileParserSmiles.h"

#include <cstdint>
#include <filesystem>

TEST_CASE("Single Aligner", "[aligner]")
{
    SECTION("Core structure is too small!")
    {
        RDKit::RWMol *mol_a = RDKit::SmilesToMol("CCCN");
        RDKit::RWMol *mol_b = RDKit::SmilesToMol("CCCO");

        coaler::SingleAligner single_aligner(5, 50);
        CHECK_THROWS_WITH(single_aligner.align_molecules_kabsch(*mol_a, *mol_b, std::nullopt),
                          "Size of core is too small!");
    };
    SECTION("Core structure is too large!")
    {
        RDKit::RWMol *mol_a = RDKit::SmilesToMol("CCCN");
        RDKit::RWMol *mol_b = RDKit::SmilesToMol("CCCO");

        coaler::SingleAligner single_aligner(1, 1);
        CHECK_THROWS_WITH(single_aligner.align_molecules_kabsch(*mol_a, *mol_b, std::nullopt),
                          "Size of core is too large!");
    };
}

