
#include <GraphMol/Substruct/SubstructMatch.h>
#include <spdlog/spdlog.h>

#include <boost/range/combine.hpp>

#include "catch2/catch.hpp"
#include "coaler/core/Forward.hpp"
#include "coaler/embedder/ConformerEmbedder.hpp"
#include "coaler/embedder/CoreSymmetryCalculator.hpp"
#include "coaler/embedder/SubstructureAnalyzer.hpp"
#include "test_helper.h"

// adapt this accoring to embedder behavior
constexpr double AVG_DEVIATION_THRESHOLD = 0;

using namespace coaler::embedder;

namespace core = coaler::core;

bool has_shared_core(const core::CoreResult& core, const RDKit::ROMOL_SPTR& mol, int confId) {
    std::vector<RDKit::MatchVectType> substructureResults;
    auto conformer = mol->getConformer(confId);
    CHECK(RDKit::SubstructMatch(*mol.get(), *core.first, substructureResults) != 0);

    RDGeom::Point3D minDiff;
    for (const auto& match : substructureResults) {
        double diffSum = 0.0;
        for (const auto& [queryId, molId] : match) {
            RDGeom::Point3D atomCoords = core.second.at(queryId);
            RDGeom::Point3D confCoords = conformer.getAtomPos(molId);
            RDGeom::Point3D diff = atomCoords - confCoords;
            diffSum += diff.length();
        }
        double avgDiff = diffSum / static_cast<double>(match.size());
        if (avgDiff <= AVG_DEVIATION_THRESHOLD) {
            return true;
        }
    }
    return false;
}

/*----------------------------------------------------------------------------------------------------------------*/

TEST_CASE("test_shared_core", "[conformer_generator_tester]") {
    auto mol1 = ROMolFromSmiles("c1ccncc1CCCO");
    auto mol2 = ROMolFromSmiles("c1c(O)cc(O)cc1O");

    RDKit::MOL_SPTR_VECT const mols = {mol1, mol2};
    auto coreResult = core::Matcher::calculateCoreMcs(mols).value();

    ConformerEmbedder embedder(coreResult.first, coreResult.second);
    embedder.embedEvenlyAcrossAllMatches(mol1, {3, 3, 3});
    // embedder.embedEvenlyAcrossAllMatches(mol2, 3);

    for (const auto& mol : mols) {
        for (int confId = 0; confId < mol->getNumConformers(); confId++) {
            CHECK(has_shared_core(coreResult, mol, confId));
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------*/
// TODO add more cases
TEST_CASE("test_ring_rotation", "[conformer_generator_tester]") {
    auto mol1 = ROMolFromSmiles("c1ccc(CC)cc1CC");
    auto mol2 = ROMolFromSmiles("c1c(CCC)cccc1");

    RDKit::MOL_SPTR_VECT const mols = {mol1, mol2};
    auto core = ROMolFromSmiles("c1ccccc1");

    RDKit::DGeomHelpers::EmbedMolecule(*core);

    RDKit::SubstructMatchParameters substructMatchParams;
    substructMatchParams.uniquify = true;
    substructMatchParams.useChirality = true;
    substructMatchParams.useQueryQueryMatches = false;
    substructMatchParams.maxMatches = 1000;
    substructMatchParams.numThreads = 1;

    auto match = RDKit::SubstructMatch(*mol1, *core, substructMatchParams).at(0);

    CoreAtomMapping molQueryCoords;
    for (int atomId = 0; atomId < core->getNumAtoms(); atomId++) {
        molQueryCoords[atomId] = core->getConformer(0).getAtomPos(atomId);
    }

    CoreAtomMapping newCoords = CoreSymmetryCalculator::getShiftedMapping(molQueryCoords, 3);
    CoreAtomMapping doubleRot = CoreSymmetryCalculator::getShiftedMapping(newCoords, 3);

    for (const auto& iter : boost::combine(molQueryCoords, doubleRot)) {
        auto oldPoint = iter.get<0>();
        auto newPoint = iter.get<1>();
        CHECK(oldPoint.second.x == newPoint.second.x);
        CHECK(oldPoint.second.y == newPoint.second.y);
        CHECK(oldPoint.second.z == newPoint.second.z);
    }
}

/*----------------------------------------------------------------------------------------------------------------*/

TEST_CASE("test_ring_symmetry_determination", "[conformer_generator_tester]") {
    std::vector<std::tuple<std::string, unsigned, unsigned>> smilesList
        = {{"C1SCNOSNC1", 1, 0}, {"c1ccccc1", 6, 1},   {"C1CCCCC1", 6, 1},  {"N1NNNNN1", 6, 1},   {"C1NCCNC1", 2, 3},
           {"C1NCNCNCN1", 4, 2}, {"C1NSOCNSO1", 2, 4}, {"C1CCCCCC1", 1, 0}, {"C1CCCCC1CCC", 1, 0}};

    for (const auto& [smiles, exp_symCount, exp_rotStep] : smilesList) {
        auto mol = *RDKit::SmilesToMol(smiles);
        auto [res_symCount, res_rotStep] = SubstructureAnalyzer::getNumberOfRingRotations(mol);
        if ((exp_symCount != res_symCount) || (exp_rotStep != res_rotStep)) {
        }
        CHECK(exp_symCount == res_symCount);
        CHECK(exp_rotStep == res_rotStep);
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
