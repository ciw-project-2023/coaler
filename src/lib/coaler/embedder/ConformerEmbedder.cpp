/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#include "ConformerEmbedder.hpp"

#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <utility>

#include "CoreSymmetryCalculator.hpp"

const unsigned seed = 42;

namespace coaler::embedder {
    ConformerEmbedder::ConformerEmbedder(RDKit::ROMol core, const CoreAtomMapping& coreMap) : m_core(std::move(core)) {
        RDKit::DGeomHelpers::EmbedParameters params;
        params.randomSeed = seed;
        params.coordMap = &coreMap;
        params.useRandomCoords = true;
        // generate Conformer with given coords for core
        RDKit::DGeomHelpers::EmbedMolecule(m_core, params);
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    bool ConformerEmbedder::embedWithFixedCore(RDKit::ROMol& mol, unsigned numConfs) {
        // match molecule and core
        std::vector<RDKit::MatchVectType> substructureResults;
        if (RDKit::SubstructMatch(mol, m_core, substructureResults) == 0) {
            return false;
        }

        // for now only use first substructure result //TODO adapt (maybe mix of different matches)
        RDKit::MatchVectType match = substructureResults.at(0);

        // determine coordinates for atoms using core conformer
        CoreAtomMapping moleculeCoreCoords;
        RDKit::Conformer coreConformer = m_core.getConformer(0);
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
        params.useRandomCoords = true;
        RDKit::DGeomHelpers::EmbedMultipleConfs(mol, numConfs, params);

        return mol.getNumConformers() == numConfs;
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    bool ConformerEmbedder::embedEvenlyAcrossSymmetryAxes(RDKit::ROMol& mol, unsigned minNofConfs,
                                                          unsigned maxNofConfs) {
        // unsigned nofSymmetryAxes = CoreSymmetryCalculator::getNofSymmetryAxes(mol);
        std::vector<RDKit::MatchVectType> substructureResults;
        if (RDKit::SubstructMatch(mol, m_core, substructureResults) == 0) {
            return false;
        }

        unsigned nofMatches = substructureResults.size();

        unsigned conformersPerMatch
            = nofMatches * minNofConfs < maxNofConfs ?(int) maxNofConfs / nofMatches : minNofConfs;
        // TODO ensure large enough

        // match all substructures
        // embedd mol with match coordinates

        return false;
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    std::vector<unsigned> ConformerEmbedder::distributeApproxEvenly(unsigned int nofMatches,
                                                                    unsigned int maxConformers) {
        std::vector<unsigned> confsForMatch(nofMatches);
        unsigned decrementPosition = maxConformers % nofMatches;
        int baseNofConfs = maxConformers / nofMatches;

        std::fill(confsForMatch.begin(),
                  confsForMatch.begin() + decrementPosition,
                  baseNofConfs + 1);

        std::fill(confsForMatch.begin() + decrementPosition + 1,
                  confsForMatch.end(),
                  (int)(maxConformers / nofMatches));

        return confsForMatch;
    }
}  // namespace coaler::embedder
