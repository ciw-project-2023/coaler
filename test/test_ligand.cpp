#include <GraphMol/SmilesParse/SmilesParse.h>

#include "catch2/catch.hpp"
#include "coaler/multialign/Forward.hpp"

using namespace coaler::multialign;

TEST_CASE("test_ligand", "[multialign]") {
    UniquePoseSet poses = {UniquePoseID(0, 0), UniquePoseID(0, 1), UniquePoseID(0, 2)};
    RDKit::RWMol mol = *RDKit::SmilesToMol("CN");
    coaler::multialign::Ligand ligand(mol, poses, 1);

    CHECK(ligand.getPoses() == poses);
    CHECK(ligand.getID() == 1);
    CHECK(ligand.getNumHeavyAtoms() == 2);
}

TEST_CASE("test_ligand_molecule", "[ligand_tester]") {
    UniquePoseSet poses = {UniquePoseID(0, 0), UniquePoseID(0, 1), UniquePoseID(0, 2)};
    RDKit::RWMol mol = *RDKit::SmilesToMol("CN");
    /*RDKit::Conformer conf(2); //TODO add conformers
    mol.addConformer(&conf);
    mol.addConformer(&conf);
    mol.addConformer(&conf);*/
    Ligand ligand(mol, poses, 1);
    CHECK(ligand.getNumPoses() == 3);
}