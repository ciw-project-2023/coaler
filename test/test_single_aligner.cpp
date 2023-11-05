#include "catch2/catch.hpp"

#include "../src/SingleAligner.hpp"

#include <cstdint>
#include <filesystem>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>

namespace test {
    using MolVec = std::vector<RDKit::RWMol*>;

    MolVec parseSmilesFile(const std::string& file_name) {
        MolVec result;

        std::string smiles;
        std::ifstream infile(file_name);
        while(std::getline(infile, smiles)) {
            RDKit::RWMol* mol =  RDKit::SmilesToMol(smiles);
            RDKit::INT_VECT conf_ids;
            RDKit::DGeomHelpers::EmbedParameters params;

            // RDKit::DGeomHelpers::EmbedMultipleConfs(*mol, conf_ids, 10, params);

            result.push_back(mol);
        }

        return result;
    }
}

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

