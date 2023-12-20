/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#include "ConformerEmbedder.hpp"

#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <spdlog/spdlog.h>

#include <boost/range/combine.hpp>
#include <utility>

#include "BucketDistributor.hpp"
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
    ConformerEmbedder::ConformerEmbedder(RDKit::ROMOL_SPTR& core, const int threads)
        : m_core(core), m_threads(threads) {}

    // TODO make shared_ptr const ref!
    void ConformerEmbedder::embedConformersWithFixedCore(const RDKit::ROMOL_SPTR& mol, unsigned numConfs) {
        // match molecule and core
        std::vector<RDKit::MatchVectType> substructureResults;
        if (RDKit::SubstructMatch(*mol, *m_core, substructureResults) == 0) {
            throw std::runtime_error("No substructure match found.");
        }

        // for now only use first substructure result //TODO adapt (maybe mix of different matches)
        const auto match = substructureResults.at(0);

        // determine coordinates for atoms using core conformer
        RDKit::Conformer coreConformer = m_core->getConformer(0);
        CoreAtomMapping moleculeCoreCoords = getAtomMappingFromMatch(match, coreConformer);
        for (const auto& matchAtom : match) {
            const int coreAtomId = matchAtom.first;
            const int molAtomId = matchAtom.second;
            const RDGeom::Point3D atomCoords = coreConformer.getAtomPos(coreAtomId);
            moleculeCoreCoords.emplace(molAtomId, atomCoords);
        }

        // embed molecule conformers
        auto params = RDKit::DGeomHelpers::srETKDGv3;
        params.randomSeed = seed;
        params.coordMap = &moleculeCoreCoords;
        params.numThreads = m_threads;
        params.useRandomCoords = true;
        RDKit::DGeomHelpers::EmbedMultipleConfs(*mol, numConfs, params);
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    bool ConformerEmbedder::embedEvenlyAcrossAllMatches(RDKit::ROMol& mol, unsigned minNofConfs, unsigned maxNofConfs) {
        // unsigned nofSymmetryAxes = CoreSymmetryCalculator::getNofSymmetryAxes(mol);
        std::vector<RDKit::MatchVectType> substructureResults;
        if (RDKit::SubstructMatch(mol, *(m_core.get()), substructureResults) == 0) {
            return false;
        }

        unsigned nofMatches = substructureResults.size();
        std::vector<unsigned> nofConformersForMatch
            = BucketDistributor::distributeApproxEvenly(nofMatches, maxNofConfs);

        if (std::any_of(nofConformersForMatch.begin(), nofConformersForMatch.end(),
                        [minNofConfs](unsigned confs) { return confs < minNofConfs; })) {
            spdlog::info(
                "Symmetry of core and/or substructure matches in structure too high for given minimum"
                "number of conformations per substructure match.");
        }
        assert(nofConformersForMatch.size() == substructureResults.size());

        for (const auto& iter : boost::combine(nofConformersForMatch, substructureResults)) {
            const unsigned nofConformers = iter.get<0>();
            const RDKit::MatchVectType& match = iter.get<1>();
            CoreAtomMapping matchCoords = getAtomMappingFromMatch(match, m_core->getConformer(0));

            auto params = RDKit::DGeomHelpers::srETKDGv3;
            params.randomSeed = seed;
            params.coordMap = &matchCoords;
            params.numThreads = m_threads;
            params.useRandomCoords = true;
            RDKit::DGeomHelpers::EmbedMultipleConfs(mol, nofConformers, params);
        }

        // TODO is this nessecary?
        return mol.getNumConformers() <= maxNofConfs;
    }

}  // namespace coaler::embedder
