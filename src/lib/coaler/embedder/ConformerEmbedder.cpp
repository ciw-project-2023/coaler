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
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <spdlog/spdlog.h>

#include <utility>

const unsigned seed = 42;
const float forceTol = 0.0135;

namespace coaler::embedder {
    ConformerEmbedder::ConformerEmbedder(coaler::core::Matcher& coreMatcher, CoreAtomMapping &coords, const int threads)
        : m_core(coreMatcher), m_threads(threads), m_coords(coords) {}

    void ConformerEmbedder::embedConformersWithFixedCore(RDKit::ROMOL_SPTR mol, unsigned numConfs) {
        spdlog::info("embedding {}", RDKit::MolToSmiles(*mol));
        spdlog::info("pattern {}", RDKit::MolToSmarts(*m_core));
        // firstMatch molecule and core
        RDKit::SubstructMatchParameters matchParams;
        auto matches = RDKit::SubstructMatch(*mol, *m_core, matchParams);
        assert(!matches.empty());

        for (auto const &match : matches) {
            CoreAtomMapping molQueryCoords;
            for (const auto &[queryId, molId] : match) {
                molQueryCoords[molId] = m_coords.at(queryId);
            }

            // embed molecule conformers
            auto params = RDKit::DGeomHelpers::srETKDGv3;;
            params.optimizerForceTol = forceTol;
            params.randomSeed = seed;
            params.coordMap = &molQueryCoords;
            params.numThreads = m_threads;

            RDKit::DGeomHelpers::EmbedMultipleConfs(*mol, numConfs, params);
            spdlog::info("embedded {} conformers.", mol->getNumConformers());

            if (mol->getNumConformers() == numConfs) {
                std::vector<std::pair<int, double>> result;
                RDKit::MMFF::MMFFOptimizeMoleculeConfs(*mol, result, m_threads);

                std::vector<unsigned> atomIds;
                for (auto const& [queryId, molId]: match) {
                    atomIds.push_back(molId);
                }

                RDKit::MolAlign::alignMolConformers(*mol, &atomIds);

                spdlog::info("optimized and aligned {} conformers.", mol->getNumConformers());

                break;
            }
        }

        assert(mol->getNumConformers() == numConfs);
    }
}  // namespace coaler::embedder
