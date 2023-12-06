
#include <GraphMol/SmilesParse/SmilesParse.h>

#include "catch2/catch.hpp"
#include "coaler/embedder/CoreSymmetryCalculator.hpp"

using namespace coaler::embedder;

TEST_CASE("test_core_symmetry", "[core_symmetry_tester]") {
    CHECK(CoreSymmetryCalculator::getNofSymmetryAxes(*RDKit::SmilesToMol("c1ccccc1")) == 12);
    CHECK(CoreSymmetryCalculator::getNofSymmetryAxes(*RDKit::SmilesToMol("c1ccncc1")) == 2);
    CHECK(CoreSymmetryCalculator::getNofSymmetryAxes(*RDKit::SmilesToMol("c1ncsc1")) == 1);
}