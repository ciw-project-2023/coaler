//
// Created by chris on 11/5/23.
//

#include "catch2/catch.hpp"
#include "../src/multialign/Forward.hpp"
#include "../src/multialign/PoseRegister.hpp"

TEST_CASE("test_add_poses_to_register", "[pose_register]") {
    coaler::multialign::PoseRegister poseRegister(1,2,3);
    coaler::multialign::PosePair pair1({0,0},{1,0});
    coaler::multialign::PosePair pair2({0,0},{2,0});
    coaler::multialign::PosePair pair3({1,0},{2,0});

    CHECK(poseRegister.addPoses(pair1, 1.7));
    CHECK(poseRegister.getHighestScoringPair() == pair1);
    CHECK(poseRegister.addPoses(pair2, 2.3));
    CHECK(poseRegister.getHighestScoringPair() == pair2);
    CHECK(poseRegister.getSize() == 2);
    CHECK(poseRegister.addPoses(pair3, 1.0));
    CHECK(poseRegister.getHighestScoringPair() == pair2);
    CHECK(poseRegister.getSize() == 3);
}

TEST_CASE("test_add_pose_to_full_register", "[pose_register]") {
    coaler::multialign::PoseRegister poseRegister(1,2,2); //limit = 2
    coaler::multialign::PosePair pair1({0,0},{1,0});
    coaler::multialign::PosePair pair2({0,0},{2,0});
    coaler::multialign::PosePair pair3({1,0},{2,0});

    CHECK(poseRegister.addPoses(pair1, 1.7));
    CHECK(poseRegister.addPoses(pair2, 2.3));
    CHECK(poseRegister.getSize() == 2);
    CHECK(!poseRegister.addPoses(pair3, 1.0));
    CHECK(poseRegister.getHighestScoringPair() == pair2);
    CHECK(poseRegister.getSize() == 2);
}