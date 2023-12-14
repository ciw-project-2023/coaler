/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#include "ConformerEmbedder.hpp"

#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <utility>

const unsigned seed = 42;

namespace coaler::embedder {
    ConformerEmbedder::ConformerEmbedder(RDKit::ROMOL_SPTR& core, const int threads)
        : m_core(core), m_threads(threads) {}

    void ConformerEmbedder::embedConformersWithFixedCore(RDKit::ROMOL_SPTR mol, unsigned numConfs) {
        // match molecule and core
        std::vector<RDKit::MatchVectType> substructureResults;
        if (RDKit::SubstructMatch(*mol, *m_core, substructureResults) == 0) {
            throw std::runtime_error("No substructure match found.");
        }

        // for now only use first substructure result //TODO adapt (maybe mix of different matches)
        const auto match = substructureResults.at(0);

        // determine coordinates for atoms using core conformer
        CoreAtomMapping moleculeCoreCoords;
        RDKit::Conformer coreConformer = m_core->getConformer(0);
        for (const auto& matchAtom : match) {
            const int coreAtomId = matchAtom.first;
            const int molAtomId = matchAtom.second;
            const RDGeom::Point3D atomCoords = coreConformer.getAtomPos(coreAtomId);
            moleculeCoreCoords.emplace(molAtomId, atomCoords);
        }

        // embed molecule conformers
        RDKit::DGeomHelpers::EmbedParameters params;
        params.randomSeed = seed;
        params.coordMap = &moleculeCoreCoords;
        params.useBasicKnowledge = true;
        params.enforceChirality = true;
        params.useSymmetryForPruning = true;
        params.useSmallRingTorsions = true;
        params.useRandomCoords = true;
        params.numThreads = m_threads;
        RDKit::DGeomHelpers::EmbedMultipleConfs(*mol, numConfs, params);
    }
}  // namespace coaler::embedder
