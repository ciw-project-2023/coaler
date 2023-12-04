
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include "catch2/catch.hpp"
#include "coaler/embedder/ConformerEmbedder.hpp"

using namespace coaler::embedder;

namespace{
    CoreAtomMapping getRandomCoreConfMapping(RDKit::ROMol& core)
    {
        RDKit::DGeomHelpers::EmbedParameters params;
        //generate Conformer with given coords for core
        RDKit::DGeomHelpers::EmbedMolecule(core, params);

        CoreAtomMapping atomMapping;

        RDKit::Conformer coreConf = core.getConformer();
        //auto coords = coreConf.

        for(int id = 0; id < core.getNumAtoms(); id++)
        {
            RDGeom::Point3D pos = coreConf.getAtomPos(id);
            atomMapping.emplace(id, pos);
        }
        return atomMapping;
    }
}

TEST_CASE("test_shared_core", "[conformer_generator_tester]") {
    RDKit::ROMol mol1 = *RDKit::SmilesToMol("c1ccccc1CCCO");
    RDKit::ROMol mol2 = *RDKit::SmilesToMol("c1c(CC)cc(CC)cc1CC");
    RDKit::ROMol core = *RDKit::SmilesToMol("c1ccccc1");

    CoreAtomMapping atomMapping = getRandomCoreConfMapping(core);

    ConformerEmbedder embedder(core, atomMapping);
    embedder.embedWithFixedCore(mol1, 10);
    embedder.embedWithFixedCore(mol2, 10);
    std::vector<RDKit::ROMol> mols = {mol1, mol2};

    for(const RDKit::ROMol& mol : mols)
    {
        std::vector<RDKit::MatchVectType> substructureResults;
        if(RDKit::SubstructMatch(mol, core , substructureResults) == 0) {
            CHECK(false);
        }
        //for now only use first substructure result //TODO adapt when changing class behavior
        RDKit::MatchVectType match = substructureResults.at(0);
        for(int id = 0; id < mol.getNumConformers(); id++)
        {
            RDKit::Conformer conf = mol.getConformer(id);

            for(const auto& matchPosition : match)
            {
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
