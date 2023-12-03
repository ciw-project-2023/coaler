//
// Created by chris on 11/20/23.
//

#pragma once
#include <unordered_map>

#include "models/Forward.hpp"

namespace coaler {
    namespace multialign {

        struct MultiAlignerResult {
          public:
            MultiAlignerResult(double score, const std::unordered_map<LigandID, PoseID>& mapping,
                               const LigandVector& ligands)
                : alignmentScore(score), poseIDsByLigandID(mapping), inputLigands(ligands) {}
            double alignmentScore;
            const std::unordered_map<LigandID, PoseID> poseIDsByLigandID;
            const std::vector<Ligand> inputLigands;
        };

    }  // namespace multialign
}  // namespace coaler