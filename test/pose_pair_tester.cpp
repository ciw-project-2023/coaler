//
// Created by follenburg on 15.11.23.
//

#include "catch2/catch.hpp"
#include "../src/multialign/Forward.hpp"
#include "../src/multialign/BasicClasses/PosePair.hpp"

TEST_CASE("test_pose_pair", "[pose_pair_tester]") {
    MultiAlign::PosePair pair1(0 ,1);
    MultiAlign::PosePair pair2(1 ,0);
    MultiAlign::PosePair pair3(2 ,1);
    MultiAlign::PosePair pair4(0 ,1);

    CHECK(pair1.getFirst() == 0);
    CHECK(pair1.getSecond() == 1);
    CHECK(pair2.getFirst() == 0);
    CHECK(pair2.getSecond() == 1);
    CHECK(pair1 == pair2);
    CHECK(!(pair1 == pair3));
    CHECK(pair1 == pair4);
    //CHECK();// check for assertion
}