/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#include "ConformerEmbedder.hpp"

#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <spdlog/spdlog.h>

const unsigned seed = 42;
const float forceTol = 0.0135;

namespace coaler::embedder {
    ConformerEmbedder::ConformerEmbedder(const core::CoreResult &result, const int threads,
                                         const bool divideConformersByMatches)
            : m_core(result), m_threads(threads), m_divideConformersByMatches(divideConformersByMatches) {}

    void ConformerEmbedder::embedConformers(const RDKit::ROMOL_SPTR &mol, unsigned numConfs) {
        // firstMatch molecule and core
        RDKit::SubstructMatchParameters substructMatchParams;
        substructMatchParams.uniquify = false;
        substructMatchParams.useChirality = false;
        substructMatchParams.useQueryQueryMatches = false;
        substructMatchParams.maxMatches = 1000;
        substructMatchParams.numThreads = m_threads;

        auto matches = RDKit::SubstructMatch(*mol, *m_core.core, substructMatchParams);
        assert(!matches.empty());

        spdlog::debug("number of Core Matches: {}", matches.size());

        unsigned matchCounter = 0;
        for (auto const &match: matches) {
            auto params = this->getEmbeddingParameters();

            std::vector<int> confs;
            if (m_divideConformersByMatches) {
                auto numConfsMatch = (matchCounter < numConfs % matches.size()) ? (int(numConfs / matches.size()) + 1)
                                                                                : (int(numConfs / matches.size()));
                confs = RDKit::DGeomHelpers::EmbedMultipleConfs(*mol, numConfsMatch, params);

                matchCounter++;
            } else {
                confs = RDKit::DGeomHelpers::EmbedMultipleConfs(*mol, numConfs, params);
            }

            std::vector<unsigned> confIds;
            for (auto const confId: confs) {
                confIds.emplace_back(confId);
            }

            std::vector<std::pair<int, double>> result;
            RDKit::MMFF::MMFFOptimizeMoleculeConfs(*mol, result, m_threads);

            // This is somewhat unintuitive:
            // We use a reference molecule (m_core.ref) to have a "real" conformer, as we cannot generate
            // chemically sensible conformers from m_core.core, as it contains smarts queries. So we get the
            // match of the query to the reference in m_core.coreToRef {'id_in_core': 'id_in_ref'}. Now we have to make a
            // list for the alignment of [(id_in_mol, id_in_ref)] because the alignment wants the atom mapping in the
            // opposite order than we get from the substruct matching (i.e. we get (queryId, molId) from substruct
            // and have to provide (molId, queryId) to the alignment.
            RDKit::MatchVectType matchReverse;
            for (const auto &[queryId, molId]: match) {
                matchReverse.emplace_back(std::make_pair(molId, m_core.coreToRef.at(queryId)));
            }

            for (auto const confId: confs) {
                auto score = RDKit::MolAlign::alignMol(*mol, *m_core.ref, confId, 0, &matchReverse);
                spdlog::debug("aligned conformer {} with score {}", confId, score);
            }
        }

        if (m_divideConformersByMatches) {
            assert(mol->getNumConformers() == numConfs);
        } else {
            assert(mol->getNumConformers() == numConfs * matches.size());
        }
    }

    RDKit::DGeomHelpers::EmbedParameters ConformerEmbedder::getEmbeddingParameters() const {
        auto params = RDKit::DGeomHelpers::srETKDGv3;
        params.randomSeed = seed;
        params.numThreads = m_threads;
        params.optimizerForceTol = forceTol;
        params.clearConfs = false;
        return params;
    }
}  // namespace coaler::embedder
