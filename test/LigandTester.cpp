//
// Created by chris on 11/5/23.
//

#include "catch2/catch.hpp"
#include "../src/multialign/Forward.hpp"
#include "../src/multialign/Ligand.hpp"

TEST_CASE("test_ligand", "[ligand_tester]") {
    std::set<MultiAlign::PoseID> poses = {2,35,23};
    MultiAlign::Ligand l(poses, 1);

    CHECK(l.getPoses() == poses);
    CHECK(l.getID() == 1);
}