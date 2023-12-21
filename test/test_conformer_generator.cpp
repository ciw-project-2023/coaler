
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <spdlog/spdlog.h>

#include <boost/range/combine.hpp>

#include "catch2/catch.hpp"
#include "coaler/core/Forward.hpp"
#include "coaler/embedder/ConformerEmbedder.hpp"
#include "coaler/embedder/CoreSymmetryCalculator.hpp"
#include "coaler/embedder/Forward.hpp"
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
        if(avgDiff <= AVG_DEVIATION_THRESHOLD) {
            return true;
        }
    }
    return false;
}

/*----------------------------------------------------------------------------------------------------------------*/

void check_confs_per_match(const RDKit::ROMOL_SPTR& mol, unsigned min, unsigned max, unsigned numMatches) {
    std::map<std::vector<double>, unsigned> counter;
    for (unsigned i = 0; i < mol->getNumConformers(); i++) {
        auto curr_conf = mol->getConformer(i);
        for (unsigned j = 0; j < curr_conf.getNumAtoms(); j++) {
            RDGeom::Point3D atomPos = curr_conf.getAtomPos(j);
            std::vector<double> atomPosVec = {atomPos.x,atomPos.y,atomPos.z};
            if (counter.find(atomPosVec) == counter.end()) {
                counter.insert({atomPosVec,1});
            }
            else {
                counter[atomPosVec]++;
            }
        }
    }
    CHECK(counter.size() == numMatches);
    unsigned map_min = counter.begin()->second;
    unsigned map_max = 0;
    for (auto entry : counter) {
        if (entry.second < map_min) {map_min = entry.second;}
        if (entry.second > map_max) {map_max = entry.second;}
    }
    // checking that max and min below/above threshold
    CHECK(map_max <= max);
    CHECK(map_min >= min);
    // checking that conformers are split equal between matches
    CHECK((map_max-map_min) <= 1);
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

/*----------------------------------------------------------------------------------------------------------------*/

TEST_CASE("test_num_conformers", "[conformer_generator_tester]") {
    auto mol1 = ROMolFromSmiles("C1CCCCC1CCCO");
    auto mol2 = ROMolFromSmiles("C1CC(C)CCC1C");


    //auto mol2 = ROMolFromSmiles("C1CCCCC1C"); this mol crashed embedder for some reason
    //auto mol2 = ROMolFromSmiles("C1CCC(C)CC1C"); this too

    //auto mol2 = ROMolFromSmiles("c1c(O)cc(O)cc1O");
    //auto mol3 = ROMolFromSmiles("c2c(CC1CC1)cc1ccccc1c2");


    //NOTE symmetric rings can only appear with murco cores (unless entire mol is just ring)
    RDKit::MOL_SPTR_VECT const mols = {mol1, mol2};
    auto coreResult = core::Matcher::calculateCoreMcs(mols).value();
    REQUIRE(coreResult.first->getNumAtoms() == 7U);
    coaler::embedder::ConformerEmbeddingParams params;
    params.maxTotalConfsPerMol = 10;
    params.maxConfsPerMatch = 10;
    params.minConfsPerMatch = 2;
    spdlog::debug(RDKit::MolToSmiles(*coreResult.first));
    ConformerEmbedder embedder(coreResult.first, coreResult.second, 1, true);
    for (const auto& mol : mols) {
        const RDKit::SubstructMatchParameters substructMatchParams = coaler::core::Matcher::getSubstructMatchParams(1);
        auto matches = RDKit::SubstructMatch(*mol, *coreResult.first, substructMatchParams);
        REQUIRE(embedder.embedEvenlyAcrossAllMatches(mol, params));
        CHECK(mol->getNumConformers() == params.maxConfsPerMatch);
        for (unsigned i = 0; i < mol->getNumConformers(); i++) {
            //check_confs_per_match(mol,params.minConfsPerMatch, params.maxConfsPerMatch, matches.size());
        }
    }

    embedder.embedEvenlyAcrossAllMatches(mol2, params);
}