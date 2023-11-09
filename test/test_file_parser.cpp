#include "catch2/catch.hpp"

#include "../src/parsers/FileParserSmiles.hpp"

#include <cstdint>
#include <filesystem>

TEST_CASE("File Parser","[parser]") {
    auto assay_targets_1806504 = FileParserSmiles::parse("test/data/AID_1806504.smi");
    REQUIRE(!assay_targets_1806504.empty());
}