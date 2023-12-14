#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/FMCS/FMCS.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>

#include <catch2/catch.hpp>
#include <coaler/io/Forward.hpp>
#include <coaler/singlealign/SingleAligner.hpp>

TEST_CASE("Singlealign", "[singlealign]") {
    auto input = coaler::io::FileParser::parse("test/data/two_mols.sdf");
    RDKit::RWMol *mol_a = input.at(0);
    RDKit::RWMol *mol_b = input.at(1);
    coaler::SingleAligner singleAligner;

    SECTION("Tanimoto Similarity") {
        SECTION("Two similar molecules") {
            double similarity = singleAligner.calculate_tanimoto_shape_similarity(*mol_a, *mol_a, -1, -1);
            CHECK(similarity == 1);
        }

        SECTION("Two different molecules") {
            double similarity = singleAligner.calculate_tanimoto_shape_similarity(*mol_a, *mol_b, -1, -1);
            // Fix for floating point precision
            similarity = static_cast<unsigned int>(similarity * 100000) / 100000.0;
            CHECK(similarity == 0.01809);
        }
    }

    SECTION("Kabsch Alignment") {
        SECTION("Two similar molecules") {
            std::vector<RDKit::ROMOL_SPTR> mols;
            mols.emplace_back(boost::make_shared<RDKit::ROMol>(*mol_a));
            mols.emplace_back(boost::make_shared<RDKit::ROMol>(*mol_a));

            RDKit::MCSResult res = RDKit::findMCS(mols);
            RDKit::ROMOL_SPTR core_structure = res.QueryMol;

            double score = singleAligner.align_molecules_kabsch(*mol_a, *mol_a, -1, -1, *core_structure);
            score = static_cast<unsigned int>(score * 100000) / 100000.0;
            CHECK(score == 0);
        }

        SECTION("Two different molecules") {
            std::vector<RDKit::ROMOL_SPTR> mols;
            mols.emplace_back(boost::make_shared<RDKit::ROMol>(*mol_a));
            mols.emplace_back(boost::make_shared<RDKit::ROMol>(*mol_b));

            RDKit::MCSResult res = RDKit::findMCS(mols);
            RDKit::ROMOL_SPTR core_structure = res.QueryMol;

            double similarity = singleAligner.align_molecules_kabsch(*mol_a, *mol_b, -1, -1, *core_structure);
            // Fix for floating point precision
            similarity = static_cast<unsigned int>(similarity * 100000) / 100000.0;
            CHECK(similarity == 1.34124);
        }

        SECTION("Core too small") {
            singleAligner = coaler::SingleAligner(20);

            std::vector<RDKit::ROMOL_SPTR> mols;
            mols.emplace_back(boost::make_shared<RDKit::ROMol>(*mol_a));
            mols.emplace_back(boost::make_shared<RDKit::ROMol>(*mol_b));

            RDKit::MCSResult res = RDKit::findMCS(mols);
            RDKit::ROMOL_SPTR core_structure = res.QueryMol;

            CHECK_THROWS(singleAligner.align_molecules_kabsch(*mol_a, *mol_b, -1, -1, *core_structure));
        }

        SECTION("Core too big") {
            singleAligner = coaler::SingleAligner(1, 20);

            std::vector<RDKit::ROMOL_SPTR> mols;
            mols.emplace_back(boost::make_shared<RDKit::ROMol>(*mol_a));
            mols.emplace_back(boost::make_shared<RDKit::ROMol>(*mol_b));

            RDKit::MCSResult res = RDKit::findMCS(mols);
            RDKit::ROMOL_SPTR core_structure = res.QueryMol;

            CHECK(singleAligner.align_molecules_kabsch(*mol_a, *mol_b, -1, -1, *core_structure));
        }
    }
}
