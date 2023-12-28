/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#pragma once
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ROMol.h>

#include "coaler/core/Forward.hpp"

#include "Forward.hpp"
#include "../multialign/Forward.hpp"

namespace coaler::embedder {

    using CoreAtomMapping = std::map<int, RDGeom::Point3D>;

    /**
     * The ConformerEmbedder class provides functionality for the generation of conformers for
     * a given molecule with contrained core coordinates.
     */
    class ConformerEmbedder {
      public:
        ConformerEmbedder(const core::CoreResult &result, const int threads, const bool divideConformersByMatches);

        /**
         * Embed an even amount of Conformers at every core match.
         * @param mol The molecule to embed.
         *
         * @note If the core has a too high symmetry, it is possible, that no embedding can be
         * performed within the given min/max constraints.
         *
         * @return True upon success.
         */
        void embedConformers(const RDKit::ROMOL_SPTR &mol, unsigned numConfs);

        //std::vector<RDKit::MatchVectType> filterMatches(const std::vector<RDKit::MatchVectType>& matches);

        static std::vector<multialign::PoseID> generateNewPosesForAssemblyLigand(const multialign::LigandPtr& ligand,
                                                           const multialign::LigandVector& targets,
                                                           const std::unordered_map<multialign::LigandID, multialign::PoseID>& conformerIDs);

        static CoreAtomMapping getLigandMcsAtomCoordsFromTargetMatch(const RDGeom::POINT3D_VECT& targetCoords,
                                                              const RDKit::MatchVectType& ligandMcsMatch,
                                                              const RDKit::MatchVectType& targetMcsMatch);

      private:
        core::CoreResult m_core;
        int m_threads;
        bool m_divideConformersByMatches;

        RDKit::DGeomHelpers::EmbedParameters getEmbeddingParameters() const;

    };
}  // namespace coaler::embedder
