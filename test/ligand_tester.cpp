//
// Created by chris on 11/5/23.
//

#include "catch2/catch.hpp"
#include "../src/multialign/Forward.hpp"
#include "../src/multialign/BasicClasses/Forward.hpp"
#include <GraphMol/SmilesParse/SmilesParse.h>

using namespace MultiAlign;
TEST_CASE("test_ligand", "[ligand_tester]") {
    UniquePoseSet poses = {
    UniquePoseIdentifier(0,0),
    UniquePoseIdentifier(0,1),
    UniquePoseIdentifier(0,2)};
    RDKit::RWMol mol = *RDKit::SmilesToMol("CN");
    MultiAlign::Ligand ligand(mol, poses, 1);

    CHECK(ligand.getPoses() == poses);
    CHECK(ligand.getID() == 1);
    CHECK(ligand.getHeavyAtomSize() == 2);
}