#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>

#include <catch2/catch.hpp>
#include <string>

#include "coaler/core/Matcher.hpp"
#include "test_helper.h"

using namespace coaler::core;

TEST_CASE("Core_constructor", "[core]") {
    auto mol1 = MolFromSmiles("c1cc(ccc1N2CCOCC2=O)N3C[C@@H](OC3=O)CNC(=O)c4ccc(s4)Cl");
    auto mol2 = MolFromSmiles("CC1(C(=O)N(C(=S)N1c2ccc(c(c2)F)C(=O)NC)c3ccc(c(c3)C(F)(F)F)C#N)C");
    RDKit::MOL_SPTR_VECT mols;
    mols.emplace_back(mol1);
    mols.emplace_back(mol2);

    auto coreScaffoldMCS = Matcher::calculateCoreMcs(mols);
    CHECK(RDKit::MolToSmarts(*coreScaffoldMCS.value().first)
          == "[#6]1:&@[#6]:&@[#6](:&@[#6]:&@[#6]:&@[#6]:&@1)-&!@[#7&R]");

    // TODO test murcko
    //
    //     auto coreScaffoldMurcko = Matcher::calculateCoreMurcko(mols);
    //     CHECK(RDKit::MolToSmarts(*coreScaffoldMurcko.value().first) == "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1");
}