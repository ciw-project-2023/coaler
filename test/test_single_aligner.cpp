#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/FMCS/FMCS.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>

#include <catch2/catch.hpp>
#include <coaler/io/Forward.hpp>
#include <coaler/singlealign/SingleAligner.hpp>

TEST_CASE("Single Aligner", "[aligner]") {
    SECTION("With H-atoms") {
        auto molecules = coaler::io::FileParser::parse("test/data/two_mols.sdf");
        RDKit::RWMol *mol_a = molecules.at(0);
        RDKit::RWMol *mol_b = molecules.at(1);

        RDKit::MolOps::addHs(*mol_a);
        RDKit::MolOps::addHs(*mol_b);

        SECTION("Use MCS") {
            spdlog::info("No core structure, start calculating MCS");
            std::vector<RDKit::ROMOL_SPTR> mols;
            mols.emplace_back(boost::make_shared<RDKit::ROMol>(*mol_a));
            mols.emplace_back(boost::make_shared<RDKit::ROMol>(*mol_b));

            RDKit::MCSResult res = RDKit::findMCS(mols);
            RDKit::ROMOL_SPTR core_structure = res.QueryMol;
            spdlog::info("MCS: " + res.SmartsString);

            SECTION("Error when core structure is too small!") {
                coaler::SingleAligner single_aligner(15, 80, true);
                CHECK_THROWS_WITH(single_aligner.align_molecules_kabsch(*mol_a, *mol_b, 0, 0, *core_structure),
                                  "Size of core is too small!");
            };
            SECTION("Warning when core structure is too large!") {
                coaler::SingleAligner single_aligner(1, 1, true);
                single_aligner.align_molecules_kabsch(*mol_a, *mol_b, 0, 0, *core_structure);
            };
            SECTION("Compare two molecules") {
                coaler::SingleAligner single_aligner(1, 80, true);
                auto result = single_aligner.align_molecules_kabsch(*mol_a, *mol_b, 0, 0, *core_structure);

                double score = static_cast<unsigned int>(result * 10000) / 10000.0;
                CHECK(score == 1.2175);
            };
        };

        SECTION("Core is only part of one molecule") {
            RDKit::RWMol *core = RDKit::SmilesToMol("OO");

            coaler::SingleAligner single_aligner(1, 80, true);

            CHECK_THROWS_WITH(single_aligner.align_molecules_kabsch(*mol_a, *mol_b, 0, 0, *core),
                              "Core is not a common core structure of molecule a and molecule b!");
        };
    };

    SECTION("Without H-atoms") {
        auto molecules = coaler::io::FileParser::parse("test/data/two_mols.sdf");
        RDKit::RWMol *mol_a = molecules.at(0);
        RDKit::RWMol *mol_b = molecules.at(1);

        SECTION("Use MCS") {
            spdlog::info("No core structure, start calculating MCS");
            std::vector<RDKit::ROMOL_SPTR> mols;
            mols.emplace_back(boost::make_shared<RDKit::ROMol>(*mol_a));
            mols.emplace_back(boost::make_shared<RDKit::ROMol>(*mol_b));

            RDKit::MCSResult res = RDKit::findMCS(mols);
            RDKit::ROMOL_SPTR core_structure = res.QueryMol;
            spdlog::info("MCS: " + res.SmartsString);

            SECTION("Error when core structure is too small!") {
                coaler::SingleAligner single_aligner(15, 80);
                CHECK_THROWS_WITH(single_aligner.align_molecules_kabsch(*mol_a, *mol_b, 0, 0, *core_structure),
                                  "Size of core is too small!");
            };
            SECTION("Warning when core structure is too large!") {
                coaler::SingleAligner single_aligner(1, 1);
                single_aligner.align_molecules_kabsch(*mol_a, *mol_b, 0, 0, *core_structure);
            };
            SECTION("Compare two molecules") {
                coaler::SingleAligner single_aligner(1, 80);
                auto result = single_aligner.align_molecules_kabsch(*mol_a, *mol_b, 0, 0, *core_structure);

                double score = static_cast<unsigned int>(result * 10000) / 10000.0;
                CHECK(score == 1.3412);
            };
        };

        SECTION("Core is only part of one molecule") {
            RDKit::RWMol *core = RDKit::SmilesToMol("OO");

            coaler::SingleAligner single_aligner(1, 80);
            CHECK_THROWS_WITH(single_aligner.align_molecules_kabsch(*mol_a, *mol_b, 0, 0, *core),
                              "Core is not a common core structure of molecule a and molecule b!");
        };
    };
}

TEST_CASE("Calc RMS", "[CalcRMS]") {
    SECTION("With H-atoms") {
        auto molecules = coaler::io::FileParser::parse("test/data/two_mols.sdf");
        RDKit::RWMol *mol_a = molecules.at(0);
        RDKit::RWMol *mol_b = molecules.at(1);

        RDKit::MolOps::addHs(*mol_a);
        RDKit::MolOps::addHs(*mol_b);

        SECTION("Use MCS") {
            spdlog::info("No core structure, start calculating MCS");
            std::vector<RDKit::ROMOL_SPTR> mols;
            mols.emplace_back(boost::make_shared<RDKit::ROMol>(*mol_a));
            mols.emplace_back(boost::make_shared<RDKit::ROMol>(*mol_b));

            RDKit::MCSResult res = RDKit::findMCS(mols);
            RDKit::ROMOL_SPTR core_structure = res.QueryMol;
            spdlog::info("MCS: " + res.SmartsString);

            SECTION("Error when core structure is too small!") {
                coaler::SingleAligner single_aligner(15, 80, true);
                CHECK_THROWS_WITH(single_aligner.calc_rms(*mol_a, *mol_b, 0, 0, *core_structure),
                                  "Size of core is too small!");
            };
            SECTION("Warning when core structure is too large!") {
                coaler::SingleAligner single_aligner(1, 1, true);
                single_aligner.calc_rms(*mol_a, *mol_b, 0, 0, *core_structure);
            };

            SECTION("Compare two molecules") {
                coaler::SingleAligner single_aligner(1, 80, true);
                auto result = single_aligner.calc_rms(*mol_a, *mol_b, 0, 0, *core_structure);

                double score = static_cast<unsigned int>(result * 10000) / 10000.0;
                // TODO: Get real score
                CHECK(score == 1.2175);
            };
        };

        SECTION("Core is only part of one molecule") {
            RDKit::RWMol *core = RDKit::SmilesToMol("OO");

            coaler::SingleAligner single_aligner(1, 80, true);

            CHECK_THROWS_WITH(single_aligner.calc_rms(*mol_a, *mol_b, 0, 0, *core),
                              "Core is not a common core structure of molecule a and molecule b!");
        };
    };
}
