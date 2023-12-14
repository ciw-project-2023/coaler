#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <spdlog/spdlog.h>

#include <catch2/catch.hpp>
#include <coaler/geometry/GeometryOptimizer.hpp>
#include <coaler/multialign/MultiAligner.hpp>
#include <coaler/singlealign/SingleAligner.hpp>
#include <iostream>

TEST_CASE("GeometryOptimizer", "[geometry]") {
    RDKit::RWMol *mol1 = RDKit::SmilesToMol("Cc1ccccc1");
    RDKit::RWMol *mol2 = RDKit::SmilesToMol("Oc1ccccc1");
    RDKit::RWMol *core = RDKit::SmilesToMol("c1ccccc1");

    RDKit::DGeomHelpers::EmbedParameters params;
    params.randomSeed = 42;
    RDKit::DGeomHelpers::EmbedMultipleConfs(*mol1, 2, params);
    RDKit::DGeomHelpers::EmbedMultipleConfs(*mol2, 2, params);
    coaler::SingleAligner singleAligner;
    std::vector<RDKit::RWMol *> mols = {mol1, mol2};
    coaler::multialign::MultiAligner aligner(mols, *core, singleAligner, 2);
    coaler::multialign::MultiAlignerResult result = aligner.alignMolecules();

    coaler::GeometryOptimizer optimizer(0.5);
    optimizer.optimize_alignment_w_icp(result);
    // coaler::multialign::MultiAlignerResult optimized_result = optimizer.get_optimized_alignment();

    auto ligands = optimizer.get_optimized_ligands();
    spdlog::info("OPTIMIZED LIGANDS {}", ligands.size());

    // TODO: remove

    const std::string file_path = "./optimized_geo.sdf";

    std::ofstream output_file(file_path);
    if (!output_file.is_open()) {
        spdlog::error("Cannot open file: {}", file_path);
        return;
    }

    boost::shared_ptr<RDKit::SDWriter> const sdf_writer(new RDKit::SDWriter(&output_file, false));
    for (const auto &entry : ligands) {
        auto conf = entry.getConformer();

        spdlog::info(std::to_string(entry.getNumConformers()));

        for (auto pos : conf.getPositions()) {
            spdlog::info("DONE {}", std::to_string(pos.x));
            spdlog::info("DONE {}", std::to_string(pos.y));
            spdlog::info("DONE {}", std::to_string(pos.z));
        }

        sdf_writer->write(entry);
    }
}
