/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#pragma once
#include <GraphMol/ROMol.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include "coaler/core/Matcher.hpp"

namespace coaler::embedder {

    using CoreAtomMapping = std::map<int, RDGeom::Point3D>;

    class ConformerEmbedder {
      public:
        // ConformerEmbedder(const RDKit::ROMol& core,
        //                   unsigned numCoreConfigs);

        explicit ConformerEmbedder(coaler::core::Matcher& coreMatcher, CoreAtomMapping& coords, int threads = 1);

        void embedConformersWithFixedCore(RDKit::ROMOL_SPTR mol, unsigned numConfs);

      private:
        coaler::core::Matcher m_core;
        CoreAtomMapping m_coords;
        int m_threads;

        std::vector<RDKit::MatchVectType> filterMatches(const std::vector<RDKit::MatchVectType>& matches);
    };
}  // namespace coaler::embedder
