/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#include "ConformerEmbedder.hpp"

#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <spdlog/spdlog.h>

#include <boost/range/combine.hpp>
#include <utility>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/MolAlign/AlignMolecules.h>

#include "BucketDistributor.hpp"
#include "CoreSymmetryCalculator.hpp"

const unsigned seed = 42;
const float forceTol = 0.0135;

/*----------------------------------------------------------------------------------------------------------------*/

namespace coaler::embedder {

    ConformerEmbedder::ConformerEmbedder(RDKit::ROMOL_SPTR &core, CoreAtomMapping coords, int threads)
            : m_core(core), m_coords(std::move(coords)), m_threads(threads) {}

    void ConformerEmbedder::embedForFirstMatch(const RDKit::ROMOL_SPTR &mol, unsigned numConfs) {
        // match molecule and core
        std::vector<RDKit::MatchVectType> substructureResults;
        if (RDKit::SubstructMatch(*mol, *m_core, substructureResults) == 0) {
            throw std::runtime_error("No substructure match found.");
        }

        // for now only use first substructure result //TODO adapt (maybe mix of different matches)
        const auto match = substructureResults.at(0);

        // determine coordinates for atoms using core conformer
        RDKit::Conformer const coreConformer = m_core->getConformer(0);

        coaler::embedder::CoreAtomMapping matchCoords;
        for (const auto &[queryId, molId]: match) {
            const RDGeom::Point3D &atomCoords = m_coords.at(queryId);
            matchCoords.emplace(molId, atomCoords);
        }

        // embed molecule conformers
        auto params = this->getEmbeddingParameters(matchCoords);
        RDKit::DGeomHelpers::EmbedMultipleConfs(*mol, numConfs, params);
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    bool
    ConformerEmbedder::embedEvenlyAcrossAllMatches(RDKit::ROMOL_SPTR &mol, unsigned minNofConfs, unsigned maxNofConfs) {
        // unsigned nofSymmetryAxes = CoreSymmetryCalculator::getNofSymmetryAxes(mol);
        std::vector<RDKit::MatchVectType> substructureResults;
        if (RDKit::SubstructMatch(*mol, *(m_core.get()), substructureResults) == 0) {
            return false;
        }

        unsigned nofMatches = substructureResults.size();
        std::vector<unsigned> nofConformersForMatch
                = BucketDistributor::distributeApproxEvenly(nofMatches, maxNofConfs);

        if (std::any_of(nofConformersForMatch.begin(), nofConformersForMatch.end(),
                        [minNofConfs](unsigned confs) { return confs < minNofConfs; })) {
            spdlog::warn(
                    "Symmetry of core and/or substructure matches in structure too high for given minimum"
                    "number of conformations per substructure match.");
        }

        assert(nofConformersForMatch.size() == substructureResults.size());

        for (const auto &iter: boost::combine(nofConformersForMatch, substructureResults)) {
            const auto nofConformers = iter.get<0>();
            const auto match = iter.get<1>();

            coaler::embedder::CoreAtomMapping matchCoords;
            for (const auto &[queryId, molId]: match) {
                const RDGeom::Point3D &atomCoords = m_coords.at(queryId);
                matchCoords.emplace(molId, atomCoords);
            }

            auto params = this->getEmbeddingParameters(matchCoords);
            RDKit::DGeomHelpers::EmbedMultipleConfs(*mol, nofConformers, params);

            spdlog::debug("embedded {} conformers into molecule: {}", nofConformers, RDKit::MolToSmiles(*mol));
        }

//        std::vector<std::pair<int, double>> result;
//        RDKit::MMFF::MMFFOptimizeMoleculeConfs(*mol, result, m_threads);

        return mol->getNumConformers() <= maxNofConfs;
    }

    RDKit::DGeomHelpers::EmbedParameters ConformerEmbedder::getEmbeddingParameters(const CoreAtomMapping& coords) {
        auto params = RDKit::DGeomHelpers::srETKDGv3;
        params.randomSeed = seed;
        params.coordMap = &coords;
        params.numThreads = m_threads;
        params.optimizerForceTol = forceTol;
        return params;
    }
}  // namespace coaler::embedder
