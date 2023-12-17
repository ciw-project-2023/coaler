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
    ConformerEmbedder::ConformerEmbedder(RDKit::ROMOL_SPTR& query, CoreAtomMapping& coords, const int threads)
        : m_core(query), m_threads(threads), m_coords(coords) {}

    void ConformerEmbedder::embedConformersWithFixedCore(RDKit::ROMOL_SPTR mol, unsigned numConfs) {
        spdlog::info("Embedding {}", RDKit::MolToSmiles(*mol));
        spdlog::info("Pattern {}", RDKit::MolToSmarts(*m_core));
        // firstMatch molecule and core
        RDKit::SubstructMatchParameters matchParams;
        auto matches = RDKit::SubstructMatch(*mol, *m_core, matchParams);
        assert(!matches.empty());

        for (auto const& match : matches) {
            CoreAtomMapping molQueryCoords;
            for (const auto& [queryId, molId] : match) {
                molQueryCoords[molId] = m_coords.at(queryId);
            }

            // embed molecule conformers
            RDKit::DGeomHelpers::EmbedParameters params;
            params = RDKit::DGeomHelpers::ETKDGv3;
            params.optimizerForceTol = forceTol;
            params.randomSeed = seed;
            params.coordMap = &molQueryCoords;
            params.useBasicKnowledge = false;
            params.enforceChirality = false;
            params.useSymmetryForPruning = false;
            params.useSmallRingTorsions = false;
            params.useRandomCoords = true;
            params.numThreads = m_threads;


            spdlog::info("OptimizerForceTol: {}", params.optimizerForceTol);



            RDKit::DGeomHelpers::EmbedMultipleConfs(*mol, numConfs, params);
            spdlog::info("Embedded {} conformers.", mol->getNumConformers());

            if (mol->getNumConformers() == numConfs) {
                std::vector<std::pair<int, double>> result;
                RDKit::MMFF::MMFFOptimizeMoleculeConfs(*mol, result, m_threads);

                spdlog::info("Optimized {} conformers.", mol->getNumConformers());

                break;
            }
        }

        assert(mol->getNumConformers() == numConfs);
    }
}  // namespace coaler::embedder
