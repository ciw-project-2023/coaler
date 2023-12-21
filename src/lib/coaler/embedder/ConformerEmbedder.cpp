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

#include "CoreSymmetryCalculator.hpp"
#include "SubstructureAnalyzer.hpp"
#include "coaler/core/Matcher.hpp"

const unsigned seed = 42;
const float forceTol = 0.0135;

namespace coaler::embedder {
    ConformerEmbedder::ConformerEmbedder(RDKit::ROMOL_SPTR &query, CoreAtomMapping &coords, const int threads,
                                         const bool divideConformersByMatches)
        : m_core(query), m_threads(threads), m_coords(coords), m_divideConformersByMatches(divideConformersByMatches) {
        for(const auto& atom : m_core->atoms())
        {
            atom->calcImplicitValence();
        }
    }

    bool ConformerEmbedder::embedEvenlyAcrossAllMatches(const RDKit::ROMOL_SPTR &mol,
                                                        const ConformerEmbeddingParams &confCountParams) {

        // firstMatch molecule and core
        RDKit::SubstructMatchParameters substructMatchParams
            = coaler::core::Matcher::getSubstructMatchParams(m_threads);

        std::vector<RDKit::MatchVectType> matches = RDKit::SubstructMatch(*mol, *m_core, substructMatchParams);
        if (matches.empty()) {
            spdlog::error("Failed to match core during conformer generation.");
            return false;
        }
        unsigned nofUniqueMatches = matches.size();
        auto [nofRotationAxes, rotationStepsize] = SubstructureAnalyzer::getNumberOfRingRotations(*m_core);
        unsigned nofTotalMatches = nofUniqueMatches * nofRotationAxes;

        // check if constraints can be upheld
        if (nofTotalMatches * confCountParams.minConfsPerMatch > confCountParams.maxTotalConfsPerMol) {
            spdlog::error("Core Symmetry or Number of Core Matches is too high for set parameters.");
            return false;  // TODO throw std::runtime_error
        }

        std::vector<unsigned> matchDistribution = CoreSymmetryCalculator::distributeApproxEvenly(
            nofTotalMatches, confCountParams.maxTotalConfsPerMol, confCountParams.maxConfsPerMatch);
        assert(matchDistribution.back() >= confCountParams.minConfsPerMatch);

        spdlog::debug("number of total core matches: {}", nofTotalMatches);

        unsigned matchPositionCounter = 0;

        // iterate over unique substructure matches
        for (auto const &match : matches) {
            // iterate over all identity core rotations
            for (int rotation = 1; rotation <= nofRotationAxes; rotation++) {
                unsigned rotationStep = rotation * rotationStepsize;
                CoreAtomMapping molQueryCoords;
                for (const auto &[queryId, molId] : match) {
                    molQueryCoords[molId] = m_coords.at(queryId);
                }

                // note: for core w/o symmetry, rotationStep becomes zero
                molQueryCoords = CoreSymmetryCalculator::getShiftedMapping(molQueryCoords, rotationStep);

                auto params = this->getEmbeddingParameters(molQueryCoords);
                params.clearConfs = false;

                unsigned nofConfsForCurrentMatch = 0;
                if (m_divideConformersByMatches) {
                    nofConfsForCurrentMatch = matchDistribution.at(matchPositionCounter);
                } else {
                    nofConfsForCurrentMatch = confCountParams.maxConfsPerMatch;
                }
                RDKit::DGeomHelpers::EmbedMultipleConfs(*mol, nofConfsForCurrentMatch, params);
            }

        }
        return true;
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    RDKit::DGeomHelpers::EmbedParameters ConformerEmbedder::getEmbeddingParameters(const CoreAtomMapping &coords) {
        auto params = RDKit::DGeomHelpers::srETKDGv3;
        params.randomSeed = seed;
        params.coordMap = &coords;
        params.numThreads = m_threads;
        params.useRandomCoords = true;
        params.optimizerForceTol = forceTol;
        return params;
    }
}  // namespace coaler::embedder
