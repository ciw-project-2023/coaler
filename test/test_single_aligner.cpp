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
    RDKit::RWMol *mol_a = RDKit::SmilesToMol("CCCN");
    RDKit::RWMol *mol_b = RDKit::SmilesToMol("CCCO");

    RDKit::ROMol molO_a = *mol_a;
    RDKit::ROMol molO_b = *mol_b;

    RDKit::DGeomHelpers::EmbedMolecule(molO_a);
    RDKit::DGeomHelpers::EmbedMolecule(molO_b);

    SECTION("Use MCS") {
        spdlog::info("No core structure, start calculating MCS");
        std::vector<RDKit::ROMOL_SPTR> mols;
        mols.emplace_back(boost::make_shared<RDKit::ROMol>(molO_a));
        mols.emplace_back(boost::make_shared<RDKit::ROMol>(molO_b));

        RDKit::MCSResult res = RDKit::findMCS(mols);
        RDKit::ROMOL_SPTR core_structure = res.QueryMol;
        spdlog::info("MCS: " + res.SmartsString);

        SECTION("Error when core structure is too small!") {
            coaler::SingleAligner single_aligner(5, 80);
            CHECK_THROWS_WITH(single_aligner.align_molecules_kabsch(molO_a, molO_b, 0, 0, *core_structure),
                              "Size of core is too small!");
        };
        SECTION("Warning when core structure is too large!") {
            coaler::SingleAligner single_aligner(1, 1);
            single_aligner.align_molecules_kabsch(molO_a, molO_b, 0, 0,
                                                  *core_structure);
        };
        SECTION("Compare two molecules") {
            coaler::SingleAligner single_aligner(1, 80);
            auto result = single_aligner.align_molecules_kabsch(molO_a, molO_b, 0, 0,
                                                                *core_structure);

            double score = static_cast<unsigned int>(result * 10000) / 10000.0;
            CHECK(score == 0.0057);
        };
    };

    SECTION("Core is only part of one molecule") {
        RDKit::RWMol *core = RDKit::SmilesToMol("CO");

        coaler::SingleAligner single_aligner(1, 80);
        CHECK_THROWS_WITH(single_aligner.align_molecules_kabsch(molO_a, molO_b, 0, 0, *core),
                          "Core is not a common core structure of molecule a and molecule b!");
    };
}

