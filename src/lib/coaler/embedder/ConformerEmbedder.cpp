/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#include "ConformerEmbedder.hpp"

#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/ForceFieldHelpers/UFF/UFF.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <spdlog/spdlog.h>

#include <utility>

const unsigned seed = 42;
const float forceTol = 0.0135;

namespace coaler::embedder {
    ConformerEmbedder::ConformerEmbedder(RDKit::ROMOL_SPTR &query, CoreAtomMapping &coords, const int threads,
                                         const bool divideConformersByMatches)
            : m_core(query), m_threads(threads), m_coords(coords),
              m_divideConformersByMatches(divideConformersByMatches) {}

    void ConformerEmbedder::embedEvenlyAcrossAllMatches(const RDKit::ROMOL_SPTR &mol, unsigned numConfs) {
        // firstMatch molecule and core
        RDKit::SubstructMatchParameters substructMatchParams;
        substructMatchParams.uniquify = false;
        substructMatchParams.useChirality = true;
        substructMatchParams.useQueryQueryMatches = false;
        substructMatchParams.maxMatches = 1000;
        substructMatchParams.numThreads = m_threads;


        auto matches = RDKit::SubstructMatch(*mol, *m_core, substructMatchParams);
        assert(!matches.empty());

        spdlog::debug("number of Core Matches: {}", matches.size());

        unsigned matchCounter = 0;
        for (auto const &match: matches) {
            CoreAtomMapping molQueryCoords;
            for (const auto &[queryId, molId]: match) {
                molQueryCoords[molId] = m_coords.at(queryId);
            }

            auto params = this->getEmbeddingParameters(molQueryCoords);

            if (m_divideConformersByMatches) {
                auto numConfsMatch = (matchCounter < numConfs % matches.size()) ? (int(numConfs / matches.size()) + 1)
                                                                                : (int(numConfs / matches.size()));
                RDKit::DGeomHelpers::EmbedMultipleConfs(*mol, numConfsMatch, params);
            } else {
                RDKit::DGeomHelpers::EmbedMultipleConfs(*mol, numConfs, params);
            }
        }

        if (m_divideConformersByMatches) {
            assert(mol->getNumConformers() == numConfs);
        } else {
            assert(mol->getNumConformers() == numConfs * matches.size());
        }
    }

    RDKit::DGeomHelpers::EmbedParameters ConformerEmbedder::getEmbeddingParameters(const CoreAtomMapping &coords) {
        auto params = RDKit::DGeomHelpers::srETKDGv3;
        params.randomSeed = seed;
        params.coordMap = &coords;
        params.numThreads = m_threads;
        params.optimizerForceTol = forceTol;
        return params;
    }
}  // namespace coaler::embedder
