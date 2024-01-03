#include <GraphMol/RDKitBase.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>

#include <catch2/catch.hpp>
#include <string>

#include "coaler/core/Matcher.hpp"
#include "coaler/io/FileParser.hpp"
#include "test_helper.h"

using namespace coaler::core;

TEST_CASE("Core_MCS", "[core]") {
    auto mol1 = MolFromSmiles("c1cc(ccc1N2CCOCC2=O)N3C[C@@H](OC3=O)CNC(=O)c4ccc(s4)Cl");
    auto mol2 = MolFromSmiles("CC1(C(=O)N(C(=S)N1c2ccc(c(c2)F)C(=O)NC)c3ccc(c(c3)C(F)(F)F)C#N)C");
    RDKit::MOL_SPTR_VECT mols;
    mols.emplace_back(mol1);
    mols.emplace_back(mol2);

    Matcher matcher(1);
    auto coreScaffoldMCS = matcher.calculateCoreMcs(mols);
    CHECK(RDKit::MolToSmarts(*coreScaffoldMCS.value().core)
          == "[#6]1:&@[#6]:&@[#6](:&@[#6]:&@[#6]:&@[#6]:&@1)-&!@[#7&R]");
}

TEST_CASE("failing_mcs", "[core]") {
    auto mol1 = MolFromSmiles("Cc1nnc2n1-c1sc3c(c1[C@@H](c1ccccc1Cl)NC2)C[C@H](C(=O)N1CCOCC1)C3");
    auto mol2 = MolFromSmiles("CCc1ccc(C2NCc3nnc(C)n3-c3sc4c(c32)CC(C(=O)N2CCOCC2)C4)cc1");
    RDKit::MOL_SPTR_VECT mols;
    mols.emplace_back(mol1);
    mols.emplace_back(mol2);

    Matcher matcher(1);
    auto coreScaffoldMCS = matcher.calculateCoreMcs(mols);
    CHECK(coreScaffoldMCS.value().core);
}

TEST_CASE("Core_Murcko", "[core]") {
    SECTION("test for a small core") {
        Matcher matcher(1);
        RDKit::MOL_SPTR_VECT smallMolPtr = coaler::io::FileParser::parse("test/data/testMurckoSmall.smi");
        auto smallCoreMurcko = matcher.calculateCoreMurcko(smallMolPtr);
        CHECK(RDKit::MolToSmarts(*smallCoreMurcko.value().core)
          == "[#6]1:&@[#6]:&@[#6]:&@[#6]:&@[#6](:&@[#6]:&@1)-&!@[#6&!R]-&!@[#6&!R](-&!@[#6&!R])-&!@[#6&!R]-&!@[#6]1:&@[#6]:&@[#6]:&@[#6]:&@[#6]:&@[#6]:&@1");
    };

    SECTION("test for a medium sized core") {
        Matcher matcher(1);
        RDKit::MOL_SPTR_VECT mediumMolPtr = coaler::io::FileParser::parse("test/data/testMurckoMedium.smi");
        auto mediumCoreMurcko = matcher.calculateCoreMurcko(mediumMolPtr);
        CHECK(RDKit::MolToSmarts(*mediumCoreMurcko.value().core)
          == "[#6]1:&@[#7,#6]:&@[#6]:&@[#6]:&@[#6]2:&@[#6]:&@1:&@[#6]:&@[#6](:&@[#6]:&@[#6]:&@2)-&!@[#6]1:&@[#7,#6]:&@[#6]:&@[#7,#6]:&@[#6]:&@[#6]:&@1");
    }

    SECTION("test for a large core") {
        Matcher matcher(1);
        RDKit::MOL_SPTR_VECT largeMolPtr = coaler::io::FileParser::parse("test/data/testMurckoLarge.smi");
        auto largeCoreMurcko = matcher.calculateCoreMurcko(largeMolPtr);
        CHECK(RDKit::MolToSmarts(*largeCoreMurcko.value().core)
          == "[#6&!R]-&!@[#6&!R](-&!@[#17,#6,#9;!R])-&!@[#8&!R]-&!@[#6&!R](-&!@[#7,#6;!R]-,=;!@[#8&!R])-&!@[#6&!R](-&!@[#6]1:&@[#6]:&@[#7,#6]:&@[#6]:&@[#7,#6]:&@[#6]:&@1)-&!@[#6&!R]-&!@[#6]1:&@[#7,#6]:&@[#6]:&@[#6]:&@[#6]2:&@[#6]:&@1:&@[#6]:&@[#6]:&@[#6]:&@[#6]:&@2");
    }
}