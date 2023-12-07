#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <coaler/io/OutputWriter.hpp>
#include <vector>

#include "catch2/catch.hpp"
using namespace coaler::multialign;
using namespace coaler::io;

TEST_CASE("Output Parser", "[io]") {
    SECTION("Add one pair to output") {
        const auto mol_a = RDKit::SmilesToMol("CCCN");
        const auto mol_b = RDKit::SmilesToMol("CCCO");

        RDKit::DGeomHelpers::EmbedMolecule(*mol_a);
        RDKit::DGeomHelpers::EmbedMolecule(*mol_b);

        const auto lig_a = Ligand(*mol_a, UniquePoseSet{0}, 0);
        const auto lig_b = Ligand(*mol_b, UniquePoseSet{0}, 0);

        auto score = 0.55;
        const MultiAlignerResult result(score,
                std::unordered_map<LigandID, PoseID>{{0, 0}, {1, 1}},
                std::vector<Ligand>{lig_a, lig_b});


        SECTION("save output in file") { OutputWriter::save_molecules_w_scores_in_file("/tmp/", result); };
    };
}
