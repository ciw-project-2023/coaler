//
// Created by follenburg on 15.11.23.
//

#include "catch2/catch.hpp"
#include "../src/multialign/Forward.hpp"
#include "../src/multialign/BasicClasses/PosePair.hpp"
#include "../src/multialign/PoseRegisterBuilder.hpp"

#include <GraphMol/SmilesParse/SmilesParse.h>


using namespace MultiAlign;

TEST_CASE("test_pose_register_builder", "[pose_register_builder_tester]") {
    PairwiseAlignment pairwiseScores;
    UniquePoseIdentifier m0p0
    Ligand l1(
            RDKit::SmilesToMol("CN"),

            
            )
    PosePair m0p0m1p0({0,0}, {1,0});
    PosePair m0p1m1p0({0,1}, {1,0});
    PosePair m0p0m1p1({0,1}, {1,1});
    pairwiseScores.emplace(m0p0m1p0, 0.8);
    pairwiseScores.emplace(m0p1m1p0, 0.5);
    pairwiseScores.emplace(m0p0m1p1, 0.1);
    PairwisePoseRegisters reg = PoseRegisterBuilder.buildPoseRegisters(
            pairwiseScores,

            )


}