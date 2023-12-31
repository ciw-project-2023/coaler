#include <GraphMol/SmilesParse/SmilesWrite.h>

#include <catch2/catch.hpp>
#include <coaler/io/FileParser.hpp>
#include <filesystem>

TEST_CASE("File Parser", "[io]") {
    SECTION("SMILES") {
        auto assay_targets_1806504 = coaler::io::FileParser::parse("test/data/AID_1806504.smi");
        REQUIRE(!assay_targets_1806504.empty());
    }
    SECTION("SDF") {
        auto molecules = coaler::io::FileParser::parse("test/data/two_mols.sdf");
        REQUIRE(molecules.size() == 2);
    }
}

TEST_CASE("Duplicates", "[io]") {
    SECTION("SMILES") {
        auto mols = coaler::io::FileParser::parse("test/data/AID_5_dup.smi");
        CHECK(mols.size() == 5);
        std::vector<std::string> smilesVec;
        for (auto mol : mols) {
            if (std::find(smilesVec.begin(), smilesVec.end(), RDKit::MolToSmiles(*mol)) == smilesVec.end()) {
                smilesVec.emplace_back(RDKit::MolToSmiles(*mol));
            } else {
                FAIL("error: duplicate not handled");
            }
        }
    }
}
