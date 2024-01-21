#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/GraphMol.h>

#include <iostream>

#include "catch2/catch.hpp"
#include "coaler/core/Forward.hpp"
#include "coaler/multialign/MultiAligner.hpp"
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

    RDKit::MOL_SPTR_VECT mols = {mol1, mol2};
    core::Matcher matcher(1);
    auto coreResult = matcher.calculateCoreMcs(mols);

    core::PairwiseMCSMap mcsMap;
    multialign::AssemblyOptimizer optimizer(mcsMap, mcsMap, 1, 0.5, 100, 1);
    multialign::MultiAligner aligner(mols, optimizer, coreResult.value(), 2);
    multialign::MultiAlignerResult result = aligner.alignMolecules();

    CHECK(result.pose_ids_by_ligand_id.size() == 2);  // ´
}
