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
        : m_core(query), m_threads(threads), m_coords(coords), m_divideConformersByMatches(divideConformersByMatches) {}

    void ConformerEmbedder::embedConformersWithFixedCore(const RDKit::ROMOL_SPTR& mol, unsigned numConfs) {
        spdlog::debug("Embedding {}", RDKit::MolToSmiles(*mol));
        spdlog::debug("Pattern {}", RDKit::MolToSmarts(*m_core));

        // firstMatch molecule and core
        RDKit::SubstructMatchParameters substructMatchParams;
        substructMatchParams.uniquify = false;
        substructMatchParams.useChirality = true;
        substructMatchParams.useQueryQueryMatches = false;
        substructMatchParams.maxMatches = 1000;
        substructMatchParams.numThreads = m_threads;

        std::vector<RDKit::MatchVectType> matches;
        unsigned nofMatches
            = RDKit::SubstructMatch(*mol, *m_core, matches, substructMatchParams.uniquify,
                                    substructMatchParams.useChirality, substructMatchParams.useQueryQueryMatches,
                                    substructMatchParams.maxMatches, substructMatchParams.numThreads);
        assert(!matches.empty());

        spdlog::debug("Number of Core Matches: {}", nofMatches);

        unsigned matchCounter = 0;
        for (auto const &match : matches) {
            CoreAtomMapping molQueryCoords;
            for (const auto &[queryId, molId] : match) {
                molQueryCoords[molId] = m_coords.at(queryId);
            }

            // embed molecule conformers
            RDKit::DGeomHelpers::EmbedParameters params;
            params = RDKit::DGeomHelpers::ETKDGv3;
            params.optimizerForceTol = forceTol;
            params.useSmallRingTorsions = true;
            params.randomSeed = seed;
            params.coordMap = &molQueryCoords;
            params.useBasicKnowledge = false;
            params.enforceChirality = false;
            params.useSymmetryForPruning = false;
            params.useSmallRingTorsions = false;
            params.useRandomCoords = true;
            params.numThreads = m_threads;
            params.clearConfs = false;

            // calculates the number of conformers for each match if m_divideConformersByMatches == true
            if (m_divideConformersByMatches) {
                unsigned nofConfsForMatch = (matchCounter < numConfs % nofMatches) ? (int(numConfs / nofMatches) + 1)
                                                                                   : (int(numConfs / nofMatches));
                RDKit::DGeomHelpers::EmbedMultipleConfs(*mol, nofConfsForMatch, params);
            } else {
                RDKit::DGeomHelpers::EmbedMultipleConfs(*mol, numConfs, params);
            }

            if ((mol->getNumConformers() == numConfs && m_divideConformersByMatches)
                || (mol->getNumConformers() == numConfs * nofMatches)) {
                std::vector<std::pair<int, double>> result;
                RDKit::MMFF::MMFFOptimizeMoleculeConfs(*mol, result, m_threads);

                spdlog::debug("Optimized {} conformers.", mol->getNumConformers());

                break;
            }
            matchCounter++;
        }

        spdlog::debug("Embedded {} conformers.", mol->getNumConformers());
        if (m_divideConformersByMatches) {
            assert(mol->getNumConformers() == numConfs);
        } else {
            spdlog::debug("mol {} numConfs {} nofMatches {}.", mol->getNumConformers(), numConfs, nofMatches);
            assert(mol->getNumConformers() == numConfs * nofMatches);
        }
    }
}  // namespace coaler::embedder
