#include "catch2/catch.hpp"
#include "coaler/multialign/models/LigandPair.hpp"

TEST_CASE("test_ligand_pair", "[multialign]") {
    coaler::multialign::LigandPair pair1(0, 1);
    coaler::multialign::LigandPair pair2(1, 0);
    coaler::multialign::LigandPair pair3(2, 1);
    coaler::multialign::LigandPair pair4(0, 1);

    CHECK(std::get<0>(pair1) == 0);
    CHECK(std::get<1>(pair1) == 1);
    CHECK(std::get<0>(pair2) == 0);
    CHECK(std::get<1>(pair2) == 1);
    CHECK(pair1 == pair2);
    CHECK(!(pair1 == pair3));
    CHECK(pair1 == pair4);
    // CHECK();// check for assertion
}