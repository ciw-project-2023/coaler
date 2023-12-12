#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <iostream>

#include <catch2/catch.hpp>
#include <coaler/geometry/GeometryOptimizer.hpp>
#include <coaler/multialign/MultiAligner.hpp>
#include <coaler/singlealign/SingleAligner.hpp>


#include <spdlog/spdlog.h>

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


    auto positions = mol2->getConformer().getPositions();
    spdlog::info("Position get");
    auto ligand = result.inputLigands.at(0).getMolecule().getConformer().getPositions();
    spdlog::info("Positions of ligand get");

    spdlog::info("Multialigner succes????");

    coaler::GeometryOptimizer optimizer(0.5);
    optimizer.optimize_alignment_w_icp(result);
}