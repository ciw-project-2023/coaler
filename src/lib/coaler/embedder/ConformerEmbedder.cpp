//
// Created by chris on 12/4/23.
//

#include "ConformerEmbedder.hpp"

#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <utility>

const unsigned seed = 42;

namespace coaler::embedder
{
    ConformerEmbedder::ConformerEmbedder(RDKit::ROMol core, const coaler::embedder::CoreAtomMapping& coreMap)
    : m_core(std::move(core)){
        RDKit::DGeomHelpers::EmbedParameters params;
        params.randomSeed = seed;
        params.coordMap = &coreMap;
        params.useRandomCoords = true;
        //generate Conformer with given coords for core
        RDKit::DGeomHelpers::EmbedMolecule(m_core, params);
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    bool ConformerEmbedder::embedWithFixedCore(const RDKit::ROMol& mol, unsigned numConfs) {
        //match molecule and core
        std::vector<RDKit::MatchVectType> substructureResults;
        if(RDKit::SubstructMatch(mol , m_core , substructureResults) == 0) {
            return false;
        }

        //for now only use first substructure result //TODO adapt (maybe mix of different matches)
        RDKit::MatchVectType match = substructureResults.at(0);

        //determine coordinates for atoms using core conformer
        CoreAtomMapping moleculeCoreCoords;
        RDKit::Conformer coreConformer = m_core.getConformer(0);
        for(const auto& matchAtom : match)
        {
            int coreAtomId = matchAtom.first;
            int molAtomId = matchAtom.second;
            RDGeom::Point3D atomCoords = coreConformer.getAtomPos(coreAtomId);
            moleculeCoreCoords.emplace(molAtomId, atomCoords);
        }
        //embedd molecule conformers
    }
}

