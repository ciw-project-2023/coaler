#include "catch2/catch.hpp"
#include "coaler/multialign/Forward.hpp"
#include "coaler/multialign/PoseRegister.hpp"

TEST_CASE("test_add_poses_to_register", "[multialign]") {
    coaler::multialign::PoseRegister poseRegister(1, 2, 2);
    coaler::multialign::PosePair pair1({0, 0}, {1, 0});
    coaler::multialign::PosePair pair2({0, 0}, {2, 0});
    coaler::multialign::PosePair pair3({1, 0}, {2, 0});

    poseRegister.addPoses(pair1, 1.7);
    CHECK(poseRegister.getSize() == 1);
    CHECK(poseRegister.getHighestScoringPair() == pair1);
    poseRegister.addPoses(pair2, 2.3);
    CHECK(poseRegister.getSize() == 2);
    CHECK(poseRegister.getHighestScoringPair() == pair2);
    poseRegister.addPoses(pair3, 1.0);
    CHECK(poseRegister.getHighestScoringPair() == pair2);
    CHECK(poseRegister.getSize() == 2);
    poseRegister.addPoses(pair3, 3.0);
    CHECK(poseRegister.getHighestScoringPair() == pair3);
    CHECK(poseRegister.getSize() == 2);
}

TEST_CASE("test_add_pose_to_full_register", "[multialign]") {
    coaler::multialign::PoseRegister poseRegister(1, 2, 2);  // limit = 2
    coaler::multialign::PosePair pair1({0, 0}, {1, 0});
    coaler::multialign::PosePair pair2({0, 0}, {2, 0});
    coaler::multialign::PosePair pair3({1, 0}, {2, 0});

    poseRegister.addPoses(pair1, 1.7);
    poseRegister.addPoses(pair2, 2.3);
    CHECK(poseRegister.getSize() == 2);
    poseRegister.addPoses(pair3, 1.0);
    CHECK(poseRegister.getHighestScoringPair() == pair2);
    CHECK(poseRegister.getSize() == 2);
}