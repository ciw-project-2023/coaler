#include <spdlog/spdlog.h>

#include "GraphMol/SmilesParse/SmartsWrite.h"
#include "GraphMol/SmilesParse/SmilesWrite.h"
#include "catch2/catch.hpp"
#include "coaler/core/Matcher.hpp"
#include "test_helper.h"

TEST_CASE("mcs contains core", "[core]") {
    auto mol1 = MolFromSmiles("c1ccccc1CC1CCCCCC1");
    auto mol2 = MolFromSmiles("c1ccccc1CCC1CCCCCC1");
    auto mol3 = MolFromSmiles("c1ccccc1C");

    auto expectedCoreMol = MolFromSmiles("c1ccccc1");

    RDKit::MOL_SPTR_VECT mols = {mol1, mol2, mol3};

    coaler::core::Matcher matcher(1);

    auto coreResult = matcher.calculateCoreMurcko(mols);
    std::string coreSmarts = RDKit::MolToSmarts(*coreResult->core);

    // ensure core is correct
    auto matches = RDKit::SubstructMatch(*expectedCoreMol, *coreResult->core);
    CHECK(!matches.empty());

    coaler::multialign::LigandVector ligands;
    coaler::multialign::LigandID id = 0;
    for (const auto& mol : mols) {
        ligands.push_back(coaler::multialign::Ligand(*mol, {}, id));
        id++;
    }
    auto pairwiseMCS = coaler::core::Matcher::calcPairwiseMCS(ligands, false, coreSmarts);

    for (const auto& mcs : pairwiseMCS) {
        std::string mcsString = get<2>(mcs.second);
        spdlog::info("shared mcs smarts: {}", mcsString);
        auto* sharedMcsMol = RDKit::SmartsToMol(mcsString);
        auto containsCore = RDKit::SubstructMatch(*sharedMcsMol, *expectedCoreMol);
        CHECK(!containsCore.empty());
    }
}
