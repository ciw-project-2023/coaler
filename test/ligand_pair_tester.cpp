//
// Created by follenburg on 15.11.23.
//

#include "catch2/catch.hpp"
// #include "coaler/multialign/Forward.hpp"
#include "coaler/multialign/models/LigandPair.hpp"

TEST_CASE("test_ligand_pair", "[ligand_pair_tester]") {
    coaler::multialign::LigandPair pair1(0, 1);
    coaler::multialign::LigandPair pair2(1, 0);
    coaler::multialign::LigandPair pair3(2, 1);
    coaler::multialign::LigandPair pair4(0, 1);

    CHECK(pair1.getFirst() == 0);
    CHECK(pair1.getSecond() == 1);
    CHECK(pair2.getFirst() == 0);
    CHECK(pair2.getSecond() == 1);
    CHECK(pair1 == pair2);
    CHECK(!(pair1 == pair3));
    CHECK(pair1 == pair4);
    // CHECK();// check for assertion
}