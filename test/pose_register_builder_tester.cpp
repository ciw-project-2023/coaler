//
// Created by follenburg on 15.11.23.
//

#include <GraphMol/SmilesParse/SmilesParse.h>

#include "catch2/catch.hpp"
#include "coaler/multialign/Forward.hpp"
#include "coaler/multialign/PoseRegister.hpp"
#include "coaler/multialign/PoseRegisterBuilder.hpp"
#include "coaler/multialign/models/Forward.hpp"

using namespace coaler::multialign;

TEST_CASE("test_pose_register_builder", "[pose_register_builder_tester]") {
    PairwiseAlignment pairwiseScores;
    UniquePoseID m0p0(0, 0);
    UniquePoseID m0p1(0, 1);
    UniquePoseID m1p0(1, 0);
    UniquePoseID m1p1(1, 1);
    Ligand l1(*RDKit::SmilesToMol("CN"), {m0p0, m0p1}, 0);
    Ligand l2(*RDKit::SmilesToMol("CO"), {m1p0, m1p1}, 1);
    LigandPair ligandPair(l1.getID(), l2.getID());
    PosePair m0p0m1p0(m0p0, m1p0);
    PosePair m0p0m1p1(m0p0, m1p1);
    PosePair m0p1m1p0(m0p1, m1p0);
    PosePair m0p1m1p1(m0p1, m1p1);
    pairwiseScores.emplace(m0p0m1p0, 0.8);
    pairwiseScores.emplace(m0p1m1p0, 0.5);
    pairwiseScores.emplace(m0p0m1p1, 0.1);
    pairwiseScores.emplace(m0p1m1p1, 0.3);

    PoseRegisterCollection collection = PoseRegisterBuilder::buildPoseRegisters(pairwiseScores, {l1, l2});

    PairwisePoseRegisters reg = collection.getAllRegisters();

    CHECK(reg.size() == 1);
    CHECK(reg.at(ligandPair)->getSize() == 4);
    CHECK(reg.at(ligandPair)->getHighestScoringPair() == m0p0m1p0);
}