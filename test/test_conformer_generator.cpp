
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <boost/range/combine.hpp>

#include "catch2/catch.hpp"
#include "coaler/embedder/ConformerEmbedder.hpp"
#include "coaler/embedder/SubstructureAnalyzer.hpp"

using namespace coaler::embedder;

namespace {
    CoreAtomMapping getRandomCoreConfMapping(RDKit::ROMol& core) {
        RDKit::DGeomHelpers::EmbedParameters params;
        // generate Conformer with given coords for core
        RDKit::DGeomHelpers::EmbedMolecule(core, params);

        CoreAtomMapping atomMapping;

        RDKit::Conformer coreConf = core.getConformer();
        // auto coords = coreConf.

        for (int id = 0; id < core.getNumAtoms(); id++) {
            RDGeom::Point3D pos = coreConf.getAtomPos(id);
            atomMapping.emplace(id, pos);
        }
        return atomMapping;
    }
}  // namespace

/*----------------------------------------------------------------------------------------------------------------*/

TEST_CASE("test_shared_core", "[conformer_generator_tester]") {
    RDKit::ROMol mol1 = *RDKit::SmilesToMol("c1ccccc1CCCO");
    RDKit::ROMol mol2 = *RDKit::SmilesToMol("c1c(CC)cc(CC)cc1CC");
    RDKit::ROMol core = *RDKit::SmilesToMol("c1ccccc1");

    CoreAtomMapping atomMapping = getRandomCoreConfMapping(core);

    ConformerEmbedder embedder(core, atomMapping);
    embedder.embedWithFixedCore(mol1, 10);
    embedder.embedWithFixedCore(mol2, 10);
    std::vector<RDKit::ROMol> mols = {mol1, mol2};

    for (const RDKit::ROMol& mol : mols) {
        std::vector<RDKit::MatchVectType> substructureResults;
        if (RDKit::SubstructMatch(mol, core, substructureResults) == 0) {
            CHECK(false);
        }
        // for now only use first substructure result //TODO adapt when changing class behavior
        RDKit::MatchVectType match = substructureResults.at(0);
        for (int id = 0; id < mol.getNumConformers(); id++) {
            RDKit::Conformer conf = mol.getConformer(id);

            for (const auto& matchPosition : match) {
                int coreAtomId = matchPosition.first;
                int molAtomId = matchPosition.second;
                RDGeom::Point3D atomCoords = core.getConformer().getAtomPos(coreAtomId);
                RDGeom::Point3D confCoords = conf.getAtomPos(molAtomId);
                RDGeom::Point3D diff = atomCoords - confCoords;
                CHECK(diff.x == 0);
                CHECK(diff.y == 0);
                CHECK(diff.z == 0);
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------*/

void check_distribution(unsigned nofMatches, unsigned maxConfs, const std::vector<unsigned>& expected_dist) {
    std::vector<unsigned> dist = ConformerEmbedder::distributeApproxEvenly(nofMatches, maxConfs);
    REQUIRE(dist.size() == expected_dist.size());
    for (const auto& itertuple : boost::combine(dist, expected_dist)) {
        CHECK(itertuple.get<0>() == itertuple.get<1>());
    }
}

/*----------------------------------------------------------------------------------------------------------------*/

TEST_CASE("validate_distribute_evenly", "[conformer_generator_tester]") {
    check_distribution(3, 7, {3, 2, 2});
    check_distribution(2, 7, {4, 3});
    check_distribution(5, 10, {2, 2, 2, 2, 2});
    check_distribution(6, 20, {4, 4, 3, 3, 3, 3});
}

/*----------------------------------------------------------------------------------------------------------------*/

TEST_CASE("test_ring_symmetry_determination", "[conformer_generator_tester]") {
    CHECK(SubstructureAnalyzer::getNumberOfRingRotations(*RDKit::SmilesToMol("c1ccccc1")) == 6);
    CHECK(SubstructureAnalyzer::getNumberOfRingRotations(*RDKit::SmilesToMol("C1CCCCC1")) == 6);
    CHECK(SubstructureAnalyzer::getNumberOfRingRotations(*RDKit::SmilesToMol("N1NNNNN1")) == 6);
    CHECK(SubstructureAnalyzer::getNumberOfRingRotations(*RDKit::SmilesToMol("C1NCCNC1")) == 2);
    CHECK(SubstructureAnalyzer::getNumberOfRingRotations(*RDKit::SmilesToMol("C1NCNCN1")) == 3);
    CHECK(SubstructureAnalyzer::getNumberOfRingRotations(*RDKit::SmilesToMol("C1NCNCNCN1")) == 4);
    CHECK(SubstructureAnalyzer::getNumberOfRingRotations(*RDKit::SmilesToMol("C1NSOCNSO1")) == 2);

    // rings of odd size are not rotation symmetric
    CHECK(SubstructureAnalyzer::getNumberOfRingRotations(*RDKit::SmilesToMol("C1CCCCCC1")) == 1);
}