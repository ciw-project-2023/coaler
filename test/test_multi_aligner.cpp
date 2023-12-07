#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <iostream>

#include "catch2/catch.hpp"
#include "coaler/multialign/Forward.hpp"
#include "coaler/multialign/MultiAligner.hpp"
#include "coaler/singlealign/SingleAligner.hpp"
#include "test_helper.h"
using namespace coaler;

TEST_CASE("basic_test", "[multialigner_tester]") {
    auto mol1 = MolFromSmiles("Cc1ccccc1");
    auto mol2 = MolFromSmiles("Oc1ccccc1");
    auto core = MolFromSmiles("c1ccccc1");

    RDKit::DGeomHelpers::EmbedParameters params;
    params.randomSeed = 42;
    RDKit::DGeomHelpers::EmbedMultipleConfs(*mol1, 2, params);
    RDKit::DGeomHelpers::EmbedMultipleConfs(*mol2, 2, params);

    SingleAligner singleAligner;
    RDKit::MOL_SPTR_VECT mols = {mol1, mol2};
    multialign::MultiAligner aligner(mols, core, singleAligner, 2);
    multialign::MultiAlignerResult result = aligner.alignMolecules();

    CHECK(result.poseIDsByLigandID.size() == 2);  // ´
}
