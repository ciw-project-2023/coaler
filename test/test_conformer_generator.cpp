
#include <GraphMol/Substruct/SubstructMatch.h>
#include <spdlog/spdlog.h>

#include "catch2/catch.hpp"
#include "coaler/core/Forward.hpp"
#include "coaler/embedder/ConformerEmbedder.hpp"
#include "test_helper.h"

//adapt this accoring to embedder behavior
constexpr double AVG_DEVIATION_THRESHOLD = 1;

using namespace coaler::embedder;

namespace core = coaler::core;

bool has_shared_core(const core::CoreResult& core, const RDKit::ROMOL_SPTR& mol, int confId){
        std::vector<RDKit::MatchVectType> substructureResults;
        auto conformer = mol->getConformer(confId);
        CHECK(RDKit::SubstructMatch(*mol.get(), *core.first, substructureResults) != 0);

    RDGeom::Point3D minDiff;
        for(const auto& match : substructureResults){
            double diffSum = 0.0;
            for (const auto& [queryId, molId] : match) {
                RDGeom::Point3D atomCoords = core.second.at(queryId);
                RDGeom::Point3D confCoords = conformer.getAtomPos(molId);
                RDGeom::Point3D diff = atomCoords - confCoords;
                diffSum += diff.length();
            }
            double avgDiff = diffSum / static_cast<double>(match.size());
            if(avgDiff <= AVG_DEVIATION_THRESHOLD) {
                return true;
            }
        }
        return false;
}


TEST_CASE("test_shared_core", "[conformer_generator_tester]") {
    auto mol1 = ROMolFromSmiles("c1ccncc1CCCO");
    auto mol2 = ROMolFromSmiles("c1c(O)cc(O)cc1O");

    RDKit::MOL_SPTR_VECT const mols = {mol1, mol2};
    auto coreResult = core::Matcher::calculateCoreMcs(mols).value();

    ConformerEmbedder embedder(coreResult.first, coreResult.second);
    embedder.embedConformersWithFixedCore(mol1, 3);
    embedder.embedConformersWithFixedCore(mol2, 3);

    for(const auto& mol : mols)
    {
        for (int confId = 0; confId <  mol->getNumConformers(); confId++){
            CHECK(has_shared_core(coreResult, mol, confId));
        }
    }
}

// TODO the embedder currently failed fot this case, I think it is because the one molecule completely contains
// the other one so the core is as big as one of the molecules. We should investigate this further.

// TEST_CASE("test_shared_core", "[conformer_generator_tester]") {
//     auto mol1 = MolFromSmiles("c1ccccc1CCCO");
//     auto mol2 = MolFromSmiles("c1c(CC)cc(CC)cc1CC");
//
//     RDKit::MOL_SPTR_VECT const mols = {mol1, mol2};
//
//     auto core = core::Matcher::calculateCoreMcs(mols);
//     ConformerEmbedder embedder(core);
//
//     auto numConfs = 10;
//     embedder.embedForFirstMatch(mol1, numConfs);
//     embedder.embedForFirstMatch(mol2, numConfs);
//
//     for (auto const& mol : mols) {
//         std::vector<RDKit::MatchVectType> substructureResults;
//         if (RDKit::SubstructMatch(*mol, *core, substructureResults) == 0) {
//             CHECK(false);
//         }
//         // for now only use first substructure result //TODO adapt when changing class behavior
//         RDKit::MatchVectType match = substructureResults.at(0);
//         for (int id = 0; id < mol->getNumConformers(); id++) {
//             RDKit::Conformer conf = mol->getConformer(id);
//
//             for (const auto& matchPosition : match) {
//                 int coreAtomId = matchPosition.first;
//                 int molAtomId = matchPosition.second;
//                 RDGeom::Point3D atomCoords = core->getConformer().getAtomPos(coreAtomId);
//                 RDGeom::Point3D confCoords = conf.getAtomPos(molAtomId);
//                 RDGeom::Point3D diff = atomCoords - confCoords;
//                 CHECK(diff.x == 0);
//                 CHECK(diff.y == 0);
//                 CHECK(diff.z == 0);
//             }
//         }
//     }
// }
