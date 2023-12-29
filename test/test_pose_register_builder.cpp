#include <GraphMol/SmilesParse/SmilesParse.h>

#include "catch2/catch.hpp"
#include "coaler/multialign/Forward.hpp"
#include "coaler/multialign/PoseRegister.hpp"
#include "coaler/multialign/PoseRegisterBuilder.hpp"
#include "coaler/multialign/models/Forward.hpp"

using namespace coaler::multialign;

TEST_CASE("test_pose_register_builder", "[multialign]") {
    PairwiseAlignments pairwiseScores;
    const UniquePoseID m0p0(0, 0);
    const UniquePoseID m0p1(0, 1);
    const UniquePoseID m1p0(1, 0);
    const UniquePoseID m1p1(1, 1);
    const Ligand l1(*RDKit::SmilesToMol("CN"), {m0p0, m0p1}, 0);
    const Ligand l2(*RDKit::SmilesToMol("CO"), {m1p0, m1p1}, 1);
    const LigandPair ligandPair(l1.getID(), l2.getID());
    const PosePair m0p0m1p0(m0p0, m1p0);
    const PosePair m0p0m1p1(m0p0, m1p1);
    const PosePair m0p1m1p0(m0p1, m1p0);
    const PosePair m0p1m1p1(m0p1, m1p1);
    pairwiseScores.emplace(m0p0m1p0, 0.8);
    pairwiseScores.emplace(m0p1m1p0, 0.5);
    pairwiseScores.emplace(m0p0m1p1, 0.1);
    pairwiseScores.emplace(m0p1m1p1, 0.3);

    SECTION("sequential") {
        const PoseRegisterCollection collection = PoseRegisterBuilder::buildPoseRegisters(pairwiseScores, {l1, l2}, 1);

        PairwisePoseRegisters reg = collection.getAllRegisters();

        CHECK(reg.size() == 1);
        CHECK(reg.at(ligandPair)->getSize()
              == Constants::POSE_REGISTER_SIZE_FACTOR * l1.getNumPoses() * l2.getNumPoses());
        CHECK(reg.at(ligandPair)->getHighestScoringPair() == m0p0m1p0);
    }
    SECTION("parallel") {
        const PoseRegisterCollection collection = PoseRegisterBuilder::buildPoseRegisters(pairwiseScores, {l1, l2}, 2);

        PairwisePoseRegisters reg = collection.getAllRegisters();

        CHECK(reg.size() == 1);
        CHECK(reg.at(ligandPair)->getSize()
              == Constants::POSE_REGISTER_SIZE_FACTOR * l1.getNumPoses() * l2.getNumPoses());
        CHECK(reg.at(ligandPair)->getHighestScoringPair() == m0p0m1p0);
    }
}