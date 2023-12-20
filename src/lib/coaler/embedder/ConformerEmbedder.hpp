/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#pragma once
#include <GraphMol/ROMol.h>
#include <GraphMol/Substruct/SubstructMatch.h>

namespace coaler::embedder {

    using CoreAtomMapping = std::map<int, RDGeom::Point3D>;

    class ConformerEmbedder {
      public:

        explicit ConformerEmbedder(RDKit::ROMOL_SPTR& query, CoreAtomMapping& coords, int threads = 1,
                                   bool divideConformersByMatches = false);

        /***
         * embedding of the molecules with their respective number of Conformers
         * @param mol
         * @param numConfs
         */
        void embedConformersWithFixedCore(RDKit::ROMOL_SPTR mol, unsigned numConfs);

      private:
        RDKit::ROMOL_SPTR m_core;
        CoreAtomMapping m_coords;
        int m_threads;
        bool m_divideConformersByMatches;

        std::vector<RDKit::MatchVectType> filterMatches(const std::vector<RDKit::MatchVectType>& matches);
    };
}  // namespace coaler::embedder
