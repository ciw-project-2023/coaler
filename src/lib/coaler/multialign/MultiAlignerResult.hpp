#pragma once
#include <unordered_map>

#include "models/Forward.hpp"

namespace coaler {
    namespace multialign {

        struct MultiAlignerResult {
            double alignment_score{};
            const std::unordered_map<LigandID, PoseID> pose_ids_by_ligand_id{};
            const std::vector<Ligand> input_ligands{};
        };

    }  // namespace multialign
}  // namespace coaler
