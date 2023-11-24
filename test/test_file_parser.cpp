#include <cstdint>
#include <filesystem>

#include "catch2/catch.hpp"
#include "io/FileParser.hpp"

TEST_CASE("File Parser", "[io]") {
    SECTION("SMILES"){
        auto assay_targets_1806504 = coaler::io::FileParser::parse("test/data/AID_1806504.smi");
        REQUIRE(!assay_targets_1806504.empty());
    }
    SECTION("SDF"){
        auto molecules = coaler::io::FileParser::parse("test/data/two_mols.sdf");
        REQUIRE(molecules.size() == 2);
    }

}
