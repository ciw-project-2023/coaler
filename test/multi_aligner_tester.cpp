//
// Created by chris on 11/9/23.
//

#include <GraphMol/DistGeomHelpers/Embedder.h>

#include <iostream>
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>

#include "../src/multialign/Forward.hpp"
#include "../src/multialign/MultiAligner.hpp"
#include "../src/singlealign/SingleAligner.hpp"
#include "catch2/catch.hpp"

using namespace coaler;

TEST_CASE("basic_test", "[multialigner_tester]") {

    RDKit::ROMol *mol1 = RDKit::SmilesToMol("Cc1ccccc1");
    RDKit::ROMol *mol2 = RDKit::SmilesToMol("Oc1ccccc1");
    RDKit::ROMol *core = RDKit::SmilesToMol("c1ccccc1");

    RDKit::DGeomHelpers::EmbedParameters params;
    params.randomSeed = 42;
    RDKit::DGeomHelpers::EmbedMultipleConfs(*mol1, 2, params);
    RDKit::DGeomHelpers::EmbedMultipleConfs(*mol2, 2, params);
    SingleAligner singleAligner;
    multialign::MultiAligner aligner({*mol1, *mol2}, *core, singleAligner, 2);
    multialign::MultiAlignerResult result = aligner.alignMolecules();
    CHECK(result.poseIDsByLigandID.size() == 2);
}