/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#include "ConformerEmbedder.hpp"

#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <spdlog/spdlog.h>

#include <boost/range/combine.hpp>
#include <utility>

#include "CoreSymmetryCalculator.hpp"

const unsigned seed = 42;

namespace {
    coaler::embedder::CoreAtomMapping getAtomMappingFromMatch(const RDKit::MatchVectType& match,
                                                              const RDKit::Conformer& matchConformer) {
        coaler::embedder::CoreAtomMapping matchCoords;
        for (const auto& matchAtom : match) {
            const int coreAtomId = matchAtom.first;
            const int molAtomId = matchAtom.second;
            const RDGeom::Point3D& atomCoords = matchConformer.getAtomPos(coreAtomId);
            matchCoords.emplace(molAtomId, atomCoords);
        }
        return matchCoords;
    }
}  // namespace

/*----------------------------------------------------------------------------------------------------------------*/

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
        CoreAtomMapping moleculeCoreCoords; //todo call getAtomMappingFromMatch()
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

    bool ConformerEmbedder::embedEvenlyAcrossAllMatches(RDKit::ROMol& mol, unsigned minNofConfs,
                                                          unsigned maxNofConfs) {
        // unsigned nofSymmetryAxes = CoreSymmetryCalculator::getNofSymmetryAxes(mol);
        std::vector<RDKit::MatchVectType> substructureResults;
        if (RDKit::SubstructMatch(mol, m_core, substructureResults) == 0) {
            return false;
        }

        unsigned nofMatches = substructureResults.size();
        std::vector<unsigned> nofConformersForMatch = distributeApproxEvenly(nofMatches, maxNofConfs);

        if(std::any_of(nofConformersForMatch.begin(), nofConformersForMatch.end(), [minNofConfs](unsigned confs){
                return confs < minNofConfs;
            })){
            spdlog::info("Symmetry of core and/or substructure matches in structure too high for given minimum"
                "number of conformations per substructure match.");
        }
        assert(nofConformersForMatch.size() == substructureResults.size());

        for(const auto& iter : boost::combine(nofConformersForMatch, substructureResults)) {
            const unsigned nofConformers = iter.get<0>();
            const RDKit::MatchVectType& match = iter.get<1>();
            CoreAtomMapping matchCoords = getAtomMappingFromMatch(match, m_core.getConformer(0));
            RDKit::DGeomHelpers::EmbedParameters params;
            params.randomSeed = seed;
            params.coordMap = &matchCoords;
            params.useRandomCoords = true;
            params.clearConfs = false;
            RDKit::DGeomHelpers::EmbedMultipleConfs(mol, nofConformers, params);
        }

        //TODO is this nessecary?
        return mol.getNumConformers() <= maxNofConfs;
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    //TODO add max per match
    std::vector<unsigned> ConformerEmbedder::distributeApproxEvenly(unsigned int nofMatches,
                                                                    unsigned int maxConformers) {
        std::vector<unsigned> confsForMatch(nofMatches);
        unsigned decrementPosition = maxConformers % nofMatches;
        unsigned baseNofConfs = maxConformers / nofMatches;

        std::fill(confsForMatch.begin(),
                  confsForMatch.begin() + decrementPosition ,
                  baseNofConfs + 1);

        std::fill(confsForMatch.begin() + decrementPosition,
                  confsForMatch.end(),
                  baseNofConfs);

        return confsForMatch;
    }
}  // namespace coaler::embedder
