/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#include "ConformerEmbedder.hpp"

#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <spdlog/spdlog.h>

const unsigned seed = 42;
const float forceTol = 0.0135;

namespace coaler::embedder {
    ConformerEmbedder::ConformerEmbedder(RDKit::ROMOL_SPTR &query, RDKit::ROMOL_SPTR &ref, const int threads,
                                         const bool divideConformersByMatches)
        : m_core(query), m_ref(ref), m_threads(threads), m_divideConformersByMatches(divideConformersByMatches) {}

    void ConformerEmbedder::embedEvenlyAcrossAllMatches(const RDKit::ROMOL_SPTR &mol, unsigned numConfs) {
        // firstMatch molecule and core
        RDKit::SubstructMatchParameters substructMatchParams;
        substructMatchParams.uniquify = false;
        substructMatchParams.useChirality = false;
        substructMatchParams.useQueryQueryMatches = false;
        substructMatchParams.maxMatches = 1000;
        substructMatchParams.numThreads = m_threads;

        auto matches = RDKit::SubstructMatch(*mol, *m_core, substructMatchParams);
        assert(!matches.empty());

        spdlog::debug("number of Core Matches: {}", matches.size());

        unsigned matchCounter = 0;
        for (auto const &match : matches) {
            RDKit::MatchVectType matchReverse;
            for (const auto &[queryId, molId] : match) {
                matchReverse.emplace_back(std::make_pair(molId, queryId));
            }

            auto params = this->getEmbeddingParameters();

            std::vector<int> confs;
            if (m_divideConformersByMatches) {
                auto numConfsMatch = (matchCounter < numConfs % matches.size()) ? (int(numConfs / matches.size()) + 1)
                                                                                : (int(numConfs / matches.size()));
                confs = RDKit::DGeomHelpers::EmbedMultipleConfs(*mol, numConfsMatch, params);
            } else {
                confs = RDKit::DGeomHelpers::EmbedMultipleConfs(*mol, numConfs, params);
            }

            std::vector<unsigned> confIds;
            for (auto const confId : confs) {
                confIds.emplace_back(confId);
            }

            std::vector<std::pair<int, double>> result;
            RDKit::MMFF::MMFFOptimizeMoleculeConfs(*mol, result, m_threads);

            for (auto const confId : confs) {
                auto score = RDKit::MolAlign::alignMol(*mol, *m_ref, confId, 0, &matchReverse);
                spdlog::debug("aligned conformer {} with score {}", confId, score);
            }
        }

        if (m_divideConformersByMatches) {
            // assert(mol->getNumConformers() == numConfs);
        } else {
            assert(mol->getNumConformers() == numConfs * matches.size());
        }
    }

    RDKit::DGeomHelpers::EmbedParameters ConformerEmbedder::getEmbeddingParameters() {
        auto params = RDKit::DGeomHelpers::srETKDGv3;
        params.randomSeed = seed;
        params.numThreads = m_threads;
        params.optimizerForceTol = forceTol;
        params.clearConfs = false;
        return params;
    }
}  // namespace coaler::embedder
