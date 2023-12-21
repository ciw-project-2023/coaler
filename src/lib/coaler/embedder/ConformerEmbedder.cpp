/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#include "ConformerEmbedder.hpp"
#include "SubstructureAnalyzer.hpp"

#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/ForceFieldHelpers/UFF/UFF.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <spdlog/spdlog.h>

#include <utility>

#include "coaler/core/Matcher.hpp"

const unsigned seed = 42;
const float forceTol = 0.0135;

namespace coaler::embedder {
    ConformerEmbedder::ConformerEmbedder(RDKit::ROMOL_SPTR &query, CoreAtomMapping &coords, const int threads,
                                         const bool divideConformersByMatches)
            : m_core(query), m_threads(threads), m_coords(coords),
              m_divideConformersByMatches(divideConformersByMatches) {}

    bool ConformerEmbedder::embedEvenlyAcrossAllMatches(const RDKit::ROMOL_SPTR &mol, const ConformerEmbeddingParams& confCountParams) {
        // firstMatch molecule and core
        const RDKit::SubstructMatchParameters substructMatchParams
            = coaler::core::Matcher::getSubstructMatchParams(m_threads);


        auto matches = RDKit::SubstructMatch(*mol, *m_core, substructMatchParams);
        assert(!matches.empty());
        unsigned nofUniqueMatches = matches.size();
        unsigned nofRotationAxes = SubstructureAnalyzer::getNumberOfRingRotations(*m_core);
        unsigned nofTotalMatches = nofUniqueMatches * nofRotationAxes;

        //check if constraints can be upheld
        if(nofTotalMatches * confCountParams.minConfsPerMatch > confCountParams.maxTotalConfsPerMol) {
            spdlog::error("Core Symmetry or Number of Core Matches is too high for set parameters.");
            return false; //TODO throw std::runtime_error
        }

        spdlog::debug("number of Core Matches: {}", nofTotalMatches);

        unsigned matchCounter = 0;
        //iterate over unique substructure matches
        for (auto const &match: matches) {

            //iterate over all identity core rotations
            for(int rotation = 1; rotation < nofRotationAxes; rotation++)
            {

            }

            CoreAtomMapping molQueryCoords;
            for (const auto &[queryId, molId]: match) {
                molQueryCoords[molId] = m_coords.at(queryId);
            }


            RDKit::DGeomHelpers::EmbedParameters embedParams = this->getEmbeddingParameters(molQueryCoords);

            if (m_divideConformersByMatches) {

                unsigned numConfsForMatch = (matchCounter < numConfs % numMatches) ? (int(numConfs / numMatches) + 1)
                                                                                   : (int(numConfs / numMatches));
                RDKit::DGeomHelpers::EmbedMultipleConfs(*mol, numConfsForMatch, embedParams);
            } else {
                RDKit::DGeomHelpers::EmbedMultipleConfs(*mol, numConfs, embedParams);
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
