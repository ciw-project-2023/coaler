//
// Created by chris on 11/5/23.
//

#include "catch2/catch.hpp"
#include "../src/multialign/Forward.hpp"
#include "../src/multialign/PoseRegister.hpp"

TEST_CASE("test_add_poses_to_register", "[pose_register]") {
    MultiAlign::PoseRegister poseRegister(1,2,2);
    MultiAlign::PosePair pair1({0,0},{1,0});
    MultiAlign::PosePair pair2({0,0},{2,0});
    MultiAlign::PosePair pair3({1,0},{2,0});

    CHECK(poseRegister.addPoses(pair1, 1.7));
    CHECK(poseRegister.getMinimumPosePairInRegister() == pair1);
    CHECK(poseRegister.addPoses(pair2, 0.3));
    CHECK(poseRegister.getMinimumPosePairInRegister() == pair2);
    CHECK(poseRegister.addPoses(pair3, 5.0));
    CHECK(poseRegister.getMinimumPosePairInRegister() == pair1);

    CHECK(poseRegister.getSize() == 2);
}