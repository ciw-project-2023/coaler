#include "catch2/catch.hpp"

#include "../src/singlealign/SingleAligner.hpp"
#include "../src/parser/FileParser.hpp"

#include <GraphMol/ROMol.h>
#include <GraphMol/FMCS/FMCS.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>

#include <cstdint>
#include <filesystem>

TEST_CASE("Single Aligner", "[aligner]") {


//    SECTION("Error when core structure is too small!") {
//        coaler::SingleAligner single_aligner(5, 80);
//        CHECK_THROWS_WITH(single_aligner.align_molecules_kabsch(molO_a, molO_b, 0, 0, std::nullopt),
//                          "Size of core is too small!");
//    };
    SECTION("Warning when core structure is too large!") {
        RDKit::RWMol *mol_a = RDKit::SmilesToMol("CCCN");
        RDKit::RWMol *mol_b = RDKit::SmilesToMol("CCCO");

        RDKit::ROMol molO_a = *mol_a;
        RDKit::ROMol molO_b = *mol_b;
        RDKit::MolOps::addHs(molO_a);
        RDKit::MolOps::addHs(molO_b);

        RDKit::DGeomHelpers::EmbedMolecule(molO_a);
        RDKit::DGeomHelpers::EmbedMolecule(molO_b);

        coaler::SingleAligner single_aligner(1, 1);
        single_aligner.align_molecules_kabsch(molO_a, molO_b, 0, 0, std::nullopt);
    };
//    SECTION("Core is only part of one molecule") {
//        RDKit::RWMol *core = RDKit::SmilesToMol("CO");
//
//        coaler::SingleAligner single_aligner(1, 80);
//        CHECK_THROWS_WITH(single_aligner.align_molecules_kabsch(molO_a, molO_b, 0, 0, *core),
//                          "Core is not a common core structure of molecule a and molecule b!");
//    };
//    SECTION("Compare two molecules") {
//        coaler::SingleAligner single_aligner(1, 80);
//        auto result = single_aligner.align_molecules_kabsch(molO_a, molO_b, 0, 0,
//                                                            std::nullopt);
//
//        double score = static_cast<unsigned int>((std::get<0>(result) * 10000)) / 10000.0;
//        CHECK(score == 0.0127);
//
//        RDKit::ROMOL_SPTR core = std::get<1>(result);
//        CHECK(RDKit::MolToSmarts(*core) == "[#6]-[#6]-[#6]");
//    };
}

