//
// Created by chris on 11/5/23.
//

#include <GraphMol/SmilesParse/SmilesParse.h>

#include "../src/multialign/Forward.hpp"
#include "catch2/catch.hpp"

using namespace coaler::multialign;
TEST_CASE("test_ligand", "[ligand_tester]") {
    UniquePoseSet poses = {UniquePoseIdentifier(0, 0), UniquePoseIdentifier(0, 1), UniquePoseIdentifier(0, 2)};
    RDKit::RWMol mol = *RDKit::SmilesToMol("CN");
    coaler::multialign::Ligand ligand(mol, poses, 1);

    CHECK(ligand.getPoses() == poses);
    CHECK(ligand.getID() == 1);
    CHECK(ligand.getHeavyAtomSize() == 2);
}

TEST_CASE("test_ligand_molecule", "[ligand_tester]") {
    UniquePoseSet poses = {UniquePoseIdentifier(0, 0), UniquePoseIdentifier(0, 1), UniquePoseIdentifier(0, 2)};
    RDKit::RWMol mol = *RDKit::SmilesToMol("CN");
    /*RDKit::Conformer conf(2); //TODO add conformers
    mol.addConformer(&conf);
    mol.addConformer(&conf);
    mol.addConformer(&conf);*/
    Ligand ligand(mol, poses, 1);
    CHECK(ligand.getNofPoses() == 3);
}