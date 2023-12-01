#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <iostream>

#include "catch2/catch.hpp"
#include "coaler/multialign/Forward.hpp"
#include "coaler/multialign/MultiAligner.hpp"
#include "coaler/singlealign/SingleAligner.hpp"

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
