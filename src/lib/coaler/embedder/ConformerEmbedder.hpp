/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#pragma once
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ROMol.h>

#include "coaler/core/Forward.hpp"
#include "coaler/multialign/Forward.hpp"

namespace coaler::embedder {
    using CoreAtomMapping = std::map<int, RDGeom::Point3D>;

    /**
     * The ConformerEmbedder class provides functionality for the generation of conformers for
     * a given molecule with contrained core coordinates.
     */
    class ConformerEmbedder {
      public:
        ConformerEmbedder(core::CoreResult result, int threads, bool divideConformersByMatches);

        /**
         * Embed an even amount of conformers at every core match.
         * @param mol The molecule to embed.
         * @param numConfs number of conformers to embed into @param mol
         *
         * @note If the core has a too high symmetry, it is possible, that no embedding can be
         * performed within the given min/max constraints.
         *
         * @return True upon success.
         */
        void embedConformers(const RDKit::ROMOL_SPTR& mol, unsigned numConfs);

        // std::vector<RDKit::MatchVectType> filterMatches(const std::vector<RDKit::MatchVectType>& matches);

        /**
         * Embed new conformers into the worst ligand of an assembly using the pairwise MCS with each target
         * @param worstLigand ligand new conformers are embedded into
         * @param targets all target ligands of the assembly
         * @param conformerIDs maps the conformerIDs to the ligands
         * @param pairwiseStrictMCSMap MCSMap of ligand pairs with strict params
         * @param pairwiseRelaxedMCSMap MCSMap of ligand pairs with relaxed params
         * @return IDs of conformers added to @param worstLigand
         */
        static std::vector<multialign::PoseID> generateNewPosesForAssemblyLigand(
            const multialign::Ligand& worstLigand, const multialign::LigandVector& targets,
            const std::unordered_map<multialign::LigandID, multialign::PoseID>& conformerIDs,
            const core::PairwiseMCSMap& pairwiseStrictMCSMap, const core::PairwiseMCSMap& pairwiseRelaxedMCSMap,
            bool enforceGeneration = false);

        /**
         * @overload
         *
         * Embed new conformers into the worst ligand of an assembly using the core structure
         * @param worstLigand ligand new conformers are embedded into
         * @param core core structure of all molecules
         * @return IDs of conformers added to @param worstLigand
         */
        static std::vector<multialign::PoseID> generateNewPosesForAssemblyLigand(const multialign::Ligand& worstLigand,
                                                                                 const core::CoreResult& core);

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
