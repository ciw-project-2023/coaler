//
// Created by follenburg on 15.11.23.
//

#include <GraphMol/SmilesParse/SmilesParse.h>

#include "../src/multialign/Forward.hpp"
#include "../src/multialign/LigandAlignmentAssembly.hpp"
#include "../src/multialign/PoseRegisterBuilder.hpp"
#include "../src/multialign/StartingAssemblyGenerator.hpp"
#include "catch2/catch.hpp"
#include <typeinfo>

using namespace coaler::multialign;

TEST_CASE("test_ligand_alignment_assembly", "[ligand_alignment_assembly_tester]") {
    UniquePoseIdentifier m0p0(0, 0);
    UniquePoseIdentifier m0p1(0, 1);
    UniquePoseIdentifier m1p0(1, 0);
    UniquePoseIdentifier m1p1(1, 1);
    Ligand l1(*RDKit::SmilesToMol("CN"), {m0p0, m0p1}, 0);
    Ligand l2(*RDKit::SmilesToMol("CO"), {m1p0, m1p1}, 1);

    PairwiseAlignment pairwiseScores;
    PosePair m0p0m1p0(m0p0, m1p0);
    PosePair m0p0m1p1(m0p0, m1p1);
    PosePair m0p1m1p0(m0p1, m1p0);
    PosePair m0p1m1p1(m0p1, m1p1);
    pairwiseScores.emplace(m0p0m1p0, 0.8);
    pairwiseScores.emplace(m0p1m1p0, 0.5);
    pairwiseScores.emplace(m0p0m1p1, 0.1);
    pairwiseScores.emplace(m0p1m1p1, 0.3);

    PoseRegisterCollection registerCollection = PoseRegisterBuilder::buildPoseRegisters
                (pairwiseScores, {l1, l2});

    LigandAlignmentAssembly assembly = StartingAssemblyGenerator::generateStartingAssembly
            (m0p0,
             registerCollection,
             {l1, l2});

    //CHECK(assembly.getPoseOfLigand(l1.getID()) == 0);
    assembly.swapPoseForLigand(0, 2);
    unsigned missingcount = assembly.getMissingLigandsCount();
    CHECK(assembly.getMissingLigandsCount() == 1);
    //CHECK(assembly.getPoseOfLigand(l1.getID()) == 2);
    assembly.incrementMissingLigandsCount();
    CHECK(assembly.getMissingLigandsCount() == missingcount+1);


    // TODO test getPoseOfLigand of of bounce for LigandID
}