#pragma once

#include <unordered_map>

#include "models/Forward.hpp"

namespace coaler {
    namespace multialign {

        struct MultiAlignerResult {
            // NOLINTBEGIN(cppcoreguidelines-pro-type-member-init, modernize-pass-by-value)
            // NOLINTNEXTLINE(readability-named-parameter)
            MultiAlignerResult(double score, const std::unordered_map<LigandID, PoseID>& mapping,
                               const LigandVector& ligands)
                : alignment_score(score), pose_ids_by_ligand_id(mapping), input_ligands(ligands) {}
            // NOLINTEND(cppcoreguidelines-pro-type-member-init, modernize-pass-by-value)
            /*--------------------------------------------------------------------------------------------------------*/

            double alignment_score;
            std::unordered_map<LigandID, PoseID> pose_ids_by_ligand_id;
            std::vector<Ligand> input_ligands;
        };

    }  // namespace multialign
}  // namespace coaler
