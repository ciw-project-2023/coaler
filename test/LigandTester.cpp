//
// Created by chris on 11/5/23.
//

#include "catch2/catch.hpp"
#include "../src/multialign/Forward.hpp"
#include "../src/multialign/Ligand.hpp"
#include <GraphMol/SmilesParse/SmilesParse.h>

TEST_CASE("test_ligand", "[ligand_tester]") {
    std::set<MultiAlign::PoseID> poses = {2,35,23};
    RDKit::RWMol mol = *RDKit::SmilesToMol("CN");
    MultiAlign::Ligand ligand(mol, poses, 1);

    CHECK(ligand.getPoses() == poses);
    CHECK(ligand.getID() == 1);
    CHECK(ligand.getHeavyAtomSize() == 2);
}