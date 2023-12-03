//
// Created by follenburg on 15.11.23.
//

#include <GraphMol/SmilesParse/SmilesParse.h>

#include "catch2/catch.hpp"
#include "coaler/multialign/Forward.hpp"
#include "coaler/multialign/LigandAlignmentAssembly.hpp"
#include "coaler/multialign/PoseRegisterBuilder.hpp"
#include "coaler/multialign/StartingAssemblyGenerator.hpp"

using namespace coaler::multialign;

TEST_CASE("test_ligand_alignment_assembly", "[ligand_alignment_assembly_tester]") {
    UniquePoseID m0p0(0, 0);
    UniquePoseID m0p1(0, 1);
    UniquePoseID m1p0(1, 0);
    UniquePoseID m1p1(1, 1);
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

    PoseRegisterCollection registerCollection = PoseRegisterBuilder::buildPoseRegisters(pairwiseScores, {l1, l2});

    LigandAlignmentAssembly assembly
        = StartingAssemblyGenerator::generateStartingAssembly(m0p0, registerCollection, {l1, l2});
    LigandID ligand1id = l1.getID();
    LigandID ligand2id = l2.getID();
    CHECK(assembly.getPoseOfLigand(ligand1id) != 0);
    CHECK(assembly.getPoseOfLigand(ligand2id) != 1);
    assembly.swapPoseForLigand(0, 2);
    unsigned missingcount = assembly.getMissingLigandsCount();
    CHECK(assembly.getMissingLigandsCount() == 1);
    CHECK((assembly.getPoseOfLigand(ligand1id) != 2));
    assembly.incrementMissingLigandsCount();
    CHECK(assembly.getMissingLigandsCount() == missingcount + 1);

    // TODO test getPoseOfLigand of of bounce for LigandID
}