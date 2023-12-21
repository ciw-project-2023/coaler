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

const unsigned seed = 42;
const float forceTol = 0.0135;

namespace{
    coaler::embedder::CoreAtomMapping getShiftedMapping(const coaler::embedder::CoreAtomMapping& map, unsigned shift) {
        coaler::embedder::CoreAtomMapping newMap;
        unsigned idMax = map.size() - 1;
        for(const auto& entry : map){
            unsigned newId = entry.first + shift;
            if(newId > idMax){
                newId = newId - idMax - 1;
            }
            newMap.emplace(newId, entry.second);
        }
        return newMap;
    }
}

namespace coaler::embedder {
    ConformerEmbedder::ConformerEmbedder(RDKit::ROMOL_SPTR &query, CoreAtomMapping &coords, const int threads,
                                         const bool divideConformersByMatches)
            : m_core(query), m_threads(threads), m_coords(coords),
              m_divideConformersByMatches(divideConformersByMatches) {}

    bool ConformerEmbedder::embedEvenlyAcrossAllMatches(const RDKit::ROMOL_SPTR &mol, const ConformerEmbeddingParams& confCountParams) {
        // firstMatch molecule and core
        RDKit::SubstructMatchParameters substructMatchParams;
        substructMatchParams.uniquify = true;
        substructMatchParams.useChirality = true;
        substructMatchParams.useQueryQueryMatches = false;
        substructMatchParams.maxMatches = 1000;
        substructMatchParams.numThreads = m_threads;


        auto matches = RDKit::SubstructMatch(*mol, *m_core, substructMatchParams);
        assert(!matches.empty());
        unsigned nofUniqueMatches = matches.size();
        auto [nofRotationAxes, rotationStepsize] = SubstructureAnalyzer::getNumberOfRingRotations(*m_core);
        unsigned nofTotalMatches = nofUniqueMatches * nofRotationAxes;

        //check if constraints can be upheld
        if(nofTotalMatches * confCountParams.minConfsPerMatch > confCountParams.maxTotalConfsPerMol) {
            spdlog::error("Core Symmetry or Number of Core Matches is too high for set parameters.");
            return false; //TODO throw std::runtime_error
        }

        spdlog::debug("number of Core Matches: {}", nofTotalMatches);
        unsigned numConfs = confCountParams.minConfsPerMatch;
        unsigned matchCounter = 0;
        //iterate over unique substructure matches
        for (auto const &match: matches) {

            //iterate over all identity core rotations
            for(int rotation = 1; rotation <= nofRotationAxes ; rotation++)
            {
                unsigned rotationStep = rotation * rotationStepsize;
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

          return true; //condition?
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
