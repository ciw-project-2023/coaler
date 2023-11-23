#include <cstdint>
#include <filesystem>

#include "../src/parser/FileParser.hpp"
#include "catch2/catch.hpp"

TEST_CASE("File Parser", "[parser]") {
    auto assay_targets_1806504 = coaler::FileParser::parse("test/data/AID_1806504.smi");
    REQUIRE(!assay_targets_1806504.empty());
}
