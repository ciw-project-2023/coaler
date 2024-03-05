#pragma once

#include <unordered_map>

#include "models/Forward.hpp"

namespace coaler {
    namespace multialign {
        /**
         * The result of a MultiAligner run.
         */
        struct MultiAlignerResult {
            /**
             * @brief Construct a new MultiAlignerResult object
             *
             * @param score The alignment score
             * @param mapping The mapping of ligand IDs to pose IDs
             * @param ligands The input ligands
             */
            // NOLINTBEGIN(cppcoreguidelines-pro-type-member-init, modernize-pass-by-value, readability-named-parameter)
            MultiAlignerResult(double score, const std::unordered_map<LigandID, PoseID>& mapping,
                               const LigandVector& ligands)
                : alignment_score(score), pose_ids_by_ligand_id(mapping), input_ligands(ligands) {}
            // NOLINTEND(cppcoreguidelines-pro-type-member-init, modernize-pass-by-value, readability-named-parameter)
            /*--------------------------------------------------------------------------------------------------------*/

            // NOLINTBEGIN(misc-non-private-member-variables-in-classes)
            double alignment_score;
            std::unordered_map<LigandID, PoseID> pose_ids_by_ligand_id;
            std::vector<Ligand> input_ligands;
            // NOLINTEND(misc-non-private-member-variables-in-classes)
        };

    }  // namespace multialign
}  // namespace coaler
