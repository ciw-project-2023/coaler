#define CATCH_CONFIG_MAIN

#include "catch2/catch.hpp"

#include "../src/SingleAligner.hpp"


#include <cstdint>
#include <GraphMol/SmilesParse/SmilesParse.h>

TEST_CASE("Aligner", "[align]") {
    RDKit::RWMol* mol_a = RDKit::SmilesToMol("CCCO");
    CHECK(1==1);
}