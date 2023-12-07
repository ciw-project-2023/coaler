#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <iostream>

#include "catch2/catch.hpp"
#include "coaler/multialign/MultiAligner.hpp"
#include "coaler/singlealign/SingleAligner.hpp"

using namespace coaler;

TEST_CASE("basic_test", "[multialign]") {
    RDKit::RWMol *mol1 = RDKit::SmilesToMol("Cc1ccccc1");
    RDKit::RWMol *mol2 = RDKit::SmilesToMol("Oc1ccccc1");
    RDKit::RWMol *core = RDKit::SmilesToMol("c1ccccc1");

    RDKit::DGeomHelpers::EmbedParameters params;
    params.randomSeed = 42;
    RDKit::DGeomHelpers::EmbedMultipleConfs(*mol1, 2, params);
    RDKit::DGeomHelpers::EmbedMultipleConfs(*mol2, 2, params);
    SingleAligner singleAligner;
    std::vector<RDKit::RWMol *> mols = {mol1, mol2};
    multialign::MultiAligner aligner(mols, *core, singleAligner, 2);
    multialign::MultiAlignerResult result = aligner.alignMolecules();
    CHECK(result.poseIDsByLigandID.size() == 2);  // ´
}
