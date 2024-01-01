#include <GraphMol/FMCS/FMCS.h>
#include <GraphMol/ROMol.h>
#include <spdlog/spdlog.h>

#include "catch2/catch.hpp"
#include "coaler/core/Matcher.hpp"
#include "test_helper.h"

TEST_CASE("test_strict_mcs", "[multialign]") {
    auto mol1 = MolFromSmiles("Cc1nnc2n1-c1sc3c(c1[C@H](c1ccccc1Cl)NC2)C[C@H](C(=O)N1CCOCC1)C3");
    auto mol2 = MolFromSmiles("Cc1nnc2n1-c1sc3c(c1[C@@H](c1ccccc1Cl)NC2)C[C@H](C(=O)N1CCOCC1)C3");
    std::vector<RDKit::ROMOL_SPTR> mols = {mol1, mol2};
    auto mcsParams = coaler::core::Matcher::getStrictMCSParams();
    auto mcsResult = RDKit::findMCS(mols, &mcsParams);
    spdlog::info("mcs: {}", mcsResult.SmartsString);
    CHECK(mcsResult.QueryMol->getNumAtoms() == mol1->getNumAtoms() - 7);

}

TEST_CASE("test_relaxed_mcs", "[multialign]") {
    auto mol1 = MolFromSmiles("Cc1nnc2n1-c1sc3c(c1[C@H](c1ccccc1Cl)NC2)C[C@H](C(=O)N1CCOCC1)C3");
    auto mol2 = MolFromSmiles("Cc1nnc2n1-c1sc3c(c1[C@@H](c1ccccc1Cl)NC2)C[C@H](C(=O)N1CCOCC1)C3");
    std::vector<RDKit::ROMOL_SPTR> mols = {mol1, mol2};
    auto mcsParams = coaler::core::Matcher::getRelaxedMCSParams();
    auto mcsResult = RDKit::findMCS(mols, &mcsParams);
    spdlog::info("mcs: {}", mcsResult.SmartsString);
    CHECK(mcsResult.QueryMol->getNumAtoms() == mol1->getNumAtoms());

}