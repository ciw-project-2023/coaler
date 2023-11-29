//
// Created by follenburg on 15.11.23.
//

#include "catch2/catch.hpp"
#include "../src/multialign/Forward.hpp"
#include "../src/multialign/LigandAlignmentAssembly.hpp"
#include "../src/multialign/PoseRegisterCollection.hpp"
#include "../src/multialign/PoseRegisterBuilder.hpp"
#include "../src/multialign/StartingAssemblyGenerator.hpp"
#include <GraphMol/SmilesParse/SmilesParse.h>

using namespace coaler::multialign;

TEST_CASE("test_ligand_alignment_assembly", "[ligand_alignment_assembly_tester]") {

    const UniquePoseIdentifier m0p0(0,0);
    const UniquePoseIdentifier m0p1(0,1);
    const UniquePoseIdentifier m1p0(1,0);
    const UniquePoseIdentifier m1p1(1,1);
    const Ligand l1(
            *RDKit::SmilesToMol("CN"),
            {m0p0, m0p1}, 0);
    const Ligand l2(
            *RDKit::SmilesToMol("CO"),
            {m1p0, m1p1}, 1);
/*
    const std::vector<Ligand> ligand{l1, l2};
    PairwiseAlignment alignmentscore;
    const PoseRegisterCollection registerCollection = PoseRegisterBuilder::buildPoseRegisters
            (alignmentscore, ligand);


    LigandAlignmentAssembly assembly = StartingAssemblyGenerator::generateStartingAssembly
            (m0p0,
             registerCollection,
             ligand);

    CHECK(assembly.getMissingLigandsCount() == 0);
    assembly.incrementMissingLigandsCount();
    CHECK(assembly.getMissingLigandsCount() == 1);
*/
    //TODO test getPoseOfLigand of of bounce for LigandID


}