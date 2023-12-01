//
// Created by follenburg on 15.11.23.
//

#include "catch2/catch.hpp"
#include "coaler/multialign/Forward.hpp"
#include "coaler/multialign/models/PosePair.hpp"

using namespace coaler::multialign;
TEST_CASE("test_pose_pair", "[pose_pair_tester]") {
    UniquePoseID pose1(0, 0);
    UniquePoseID pose2(1, 0);
    UniquePoseID pose3(2, 0);

    PosePair pair1(pose1, pose2);
    PosePair pair2(pose2, pose1);
    PosePair pair3(pose1, pose3);
    PosePair pair4(pose2, pose3);

    CHECK(pair1.getFirst() == pose1);
    CHECK(pair1.getSecond() == pose2);
    CHECK(pair2.getFirst() == pose1);
    CHECK(pair2.getSecond() == pose2);
    CHECK(pair1 == pair2);
    CHECK(!(pair1 == pair3));
    CHECK(!(pair4 == pair3));
    // CHECK();// check for assertion
}