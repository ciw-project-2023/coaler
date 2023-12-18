//
// Created by malte on 11/27/23.
//

#include "catch2/catch.hpp"
#include "coaler/multialign/Forward.hpp"
#include "coaler/multialign/models/UniquePoseID.hpp"

using namespace coaler::multialign;

TEST_CASE("test_getter_functions", "[unique_pose_identifier_tester]") {
    UniquePoseID pose1(0, 0);
    CHECK(pose1.getLigandId() == 0);
    CHECK(pose1.getLigandInternalPoseId() == 0);

    UniquePoseID pose2(1, 1);
    CHECK(pose2.getLigandId() == 1);
    CHECK(pose2.getLigandInternalPoseId() == 1);
}

TEST_CASE("test_operator_functions", "[unique_pose_identifier_tester]") {
    UniquePoseID pose1(0, 0);
    UniquePoseID pose2(0, 0);
    UniquePoseID pose3(0, 1);
    UniquePoseID pose4(1, 0);
    UniquePoseID pose5(1, 1);
    CHECK(pose1 == pose2);
    CHECK(pose1 != pose3);
    CHECK(pose1 != pose4);
    CHECK(pose1 != pose5);
    CHECK(pose3 != pose4);
    CHECK(pose3 != pose5);
    CHECK(pose4 != pose5);
    CHECK_FALSE(pose1 == pose3);
    CHECK_FALSE(pose1 != pose1);

    CHECK((pose1 == pose2) == (pose2 == pose1));
    CHECK((pose3 == pose4) == (pose4 == pose3));
    CHECK_FALSE((pose3 != pose4) == (pose4 == pose3));
    CHECK_FALSE((pose3 == pose4) == (pose4 != pose3));

    CHECK(pose1 < pose3);
    CHECK(pose1 < pose4);
    CHECK(pose1 < pose5);
    CHECK(pose3 < pose4);
    CHECK(pose3 < pose5);
    CHECK(pose4 < pose5);
    CHECK_FALSE(pose1 < pose1);
    CHECK_FALSE(pose3 < pose1);

    CHECK(pose3 > pose1);
    CHECK(pose4 > pose1);
    CHECK(pose5 > pose1);
    CHECK(pose4 > pose3);
    CHECK(pose5 > pose3);
    CHECK(pose5 > pose4);
    CHECK_FALSE(pose1 > pose1);
    CHECK_FALSE(pose1 > pose3);
}

TEST_CASE("test_string_function", "[unique_pose_identifier_tester]") {
    UniquePoseID pose1(0, 0);
    UniquePoseID pose2(1, 1);

    CHECK(pose1.toString() == "0-0");
    CHECK(pose2.toString() == "1-1");
}

TEST_CASE("test_hash_operator", "[unique_pose_identifier_tester]") {
    UniquePoseID pose1(0, 0);
    UniquePoseID pose2(1, 0);
    UniquePoseID pose3(1, 1);

    std::size_t hashed_pose_1 = std::hash<std::string>{}("0-0");
    std::size_t hashed_pose_2 = std::hash<std::string>{}("1-0");
    std::size_t hashed_pose_none = std::hash<std::string>{}("0-1");

    CHECK(hashed_pose_1 == UniquePoseIdentifierHash()(pose1));
    CHECK_FALSE(hashed_pose_1 == UniquePoseIdentifierHash()(pose2));
    CHECK_FALSE(hashed_pose_2 == UniquePoseIdentifierHash()(pose3));
    CHECK_FALSE(hashed_pose_none == UniquePoseIdentifierHash()(pose2));
}