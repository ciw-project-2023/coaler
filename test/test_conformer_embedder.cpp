
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include "GraphMol/SmilesParse/SmilesWrite.h"
#include "catch2/catch.hpp"
#include "coaler/core/Forward.hpp"
#include "coaler/embedder/ConformerEmbedder.hpp"
#include "coaler/embedder/SubstructureAnalyzer.hpp"
#include "test_helper.h"

using namespace coaler::embedder;
namespace core = coaler::core;

TEST_CASE("test_mcs", "[conformer_generator_tester]") {
    SECTION("test with one match") {
        auto mol1 = ROMolFromSmiles("c1ccncc1CCCO");
        auto mol2 = ROMolFromSmiles("c1c(O)cc(O)cc1O");

        core::Matcher matcher(1);

        RDKit::MOL_SPTR_VECT mols = {mol1, mol2};
        auto core = matcher.calculateCoreMcs(mols).value();

        ConformerEmbedder embedder(core, 1, true);
        embedder.embedConformers(mol1, 10);

        CHECK(RDKit::MolToSmarts(*core.core) == "[#6]1:&@[#6](-&!@[#8,#6;!R]):&@[#6]:&@[#6,#7]:&@[#6]:&@[#6]:&@1");
        CHECK(mol1->getNumConformers() == 10);
    }

    SECTION("with 7 matches and enough conformers") {
        auto mol1 = ROMolFromSmiles("CC1=C2C(C(C)=CC3=C2C4=C(CC(CCC)C4)C=C3C)=CC=C1");
        auto mol2 = ROMolFromSmiles("C1(CC(C=CC2)=C2C3)=C3C=CC=C1");

        core::Matcher matcher(1);

        RDKit::MOL_SPTR_VECT mols = {mol1, mol2};
        auto core = matcher.calculateCoreMcs(mols).value();

        ConformerEmbedder embedder(core, 1, true);
        embedder.embedConformers(mol1, 70);

        CHECK(RDKit::MolToSmarts(*core.core)
              == "[#6]12-,:;@[#6]-,:;@[#6]=,:;@[#6]-,:;@[#6]-,:;@[#6]:&@1:&@[#6]:&@[#6]:&@[#6]:&@[#6]:&@2");
        CHECK(mol1->getNumConformers() == 70);
    }

    SECTION("with 7 matches and to little conformers") {
        auto mol1 = ROMolFromSmiles("CC1=C2C(C(C)=CC3=C2C4=C(CC(CCC)C4)C=C3C)=CC=C1");
        auto mol2 = ROMolFromSmiles("C1(CC(C=CC2)=C2C3)=C3C=CC=C1");

        core::Matcher matcher(1);

        RDKit::MOL_SPTR_VECT mols = {mol1, mol2};
        auto core = matcher.calculateCoreMcs(mols).value();

        ConformerEmbedder embedder(core, 1, true);
        embedder.embedConformers(mol1, 30);

        CHECK(RDKit::MolToSmarts(*core.core)
              == "[#6]12-,:;@[#6]-,:;@[#6]=,:;@[#6]-,:;@[#6]-,:;@[#6]:&@1:&@[#6]:&@[#6]:&@[#6]:&@[#6]:&@2");
        CHECK(mol1->getNumConformers() == 40);
    }
}

TEST_CASE("test_mcs_match_mapping", "[conformer_generator_tester]") {
    RDGeom::POINT3D_VECT points;
    for (int i = 0; i <= 5; i++) {
        points.emplace_back(i, 0, 0);
    }
    const RDKit::MatchVectType targetMatch = {{0, 2}, {1, 3}, {2, 4}};
    const RDKit::MatchVectType ligandMatch = {{0, 2}, {1, 3}, {2, 4}};
    auto mol1 = ROMolFromSmiles("CC1=C2C(C(C)=CC3=C2C4=C(CC(CCC)C4)C=C3C)=CC=C1");
    auto mol2 = ROMolFromSmiles("C1(CC(C=CC2)=C2C3)=C3C=CC=C1");

    core::Matcher matcher(1);

    RDKit::MOL_SPTR_VECT mols = {mol1, mol2};
    auto core = matcher.calculateCoreMcs(mols).value();

    ConformerEmbedder embedder(core, 1, true);
    embedder.embedConformers(mol1, 30);
    CoreAtomMapping ligandCoords
        = embedder.getLigandMcsAtomCoordsFromTargetMatch(points, ligandMatch, targetMatch);

    CHECK(ligandCoords.count(2) == 1);
    CHECK(ligandCoords.count(3) == 1);
    CHECK(ligandCoords.count(4) == 1);

    CHECK(ligandCoords.at(2).x == 2);
    CHECK(ligandCoords.at(3).x == 3);
    CHECK(ligandCoords.at(4).x == 4);
}
