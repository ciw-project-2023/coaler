/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#pragma once
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ROMol.h>

#include "coaler/core/Forward.hpp"
#include "coaler/multialign/models/Forward.hpp"
/**
 * @file ConformerEmbedder.hpp
 * @brief This file contains the ConformerEmbedder class which is used to generate conformers for a given molecule
 */
namespace coaler::embedder {
    using CoreAtomMapping = std::map<int, RDGeom::Point3D>;

    /**
     * The ConformerEmbedder class provides functionality for the generation of conformers for
     * a given molecule with contrained core coordinates.
     */
    class ConformerEmbedder {
      public:
        ConformerEmbedder(core::CoreResult result, int threads, bool divideConformersByMatches = true);

        /**
         * Embeds conformers into a molecule
         * @param mol molecule to embed conformers into
         * @param numConfs number of conformers to embed
         */
        void embedConformers(const RDKit::ROMOL_SPTR& mol, unsigned numConfs);

        // std::vector<RDKit::MatchVectType> filterMatches(const std::vector<RDKit::MatchVectType>& matches);

        /**
         * Embeds conformers into a molecule
         *
         * @param mol molecule to embed conformers into
         * @param numConfs number of conformers to embed
         * @param matches matches to use for embedding
         *
         * @return vector of new poses
         */
        static std::vector<multialign::PoseID> generateNewPosesForAssemblyLigand(
            const multialign::Ligand& worstLigand, const multialign::LigandVector& targets,
            const std::unordered_map<multialign::LigandID, multialign::PoseID>& conformerIDs,
            const core::PairwiseMCSMap& pairwiseStrictMCSMap, const core::PairwiseMCSMap& pairwiseRelaxedMCSMap,
            bool enforceGeneration = false);

        /**
         * @overload generateNewPosesForAssemblyLigand
         *
         * @param worstLigand
         * @param numConfs
         *
         * @return vector of new poses
         */
        std::vector<multialign::PoseID> generateNewPosesForAssemblyLigand(const multialign::Ligand& worstLigand,
                                                                          const unsigned numConfs);

        /**
         * Get the coordinates of the ligand MCS atoms from the target match
         *
         * @param targetCoords
         * @param ligandMcsMatch
         * @param targetMcsMatch
         * @return CoreAtomMapping
         */
        static CoreAtomMapping getLigandMcsAtomCoordsFromTargetMatch(const RDGeom::POINT3D_VECT& targetCoords,
                                                                     const RDKit::MatchVectType& ligandMcsMatch,
                                                                     const RDKit::MatchVectType& targetMcsMatch);

      private:
        core::CoreResult m_core;
        int m_threads;
        bool m_divideConformersByMatches;

        [[nodiscard]] RDKit::DGeomHelpers::EmbedParameters getEmbeddingParameters() const;
    };
}  // namespace coaler::embedder
