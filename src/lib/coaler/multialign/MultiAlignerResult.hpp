/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#pragma once
#include <unordered_map>
#include <vector>

#include "Forward.hpp"

namespace coaler::multialign {

    struct MultiAlignerResult {
      public:
        MultiAlignerResult(double score, const std::unordered_map<LigandID, PoseID>& mapping,
                           const LigandVector& ligands)
            : alignmentScore(score), poseIDsByLigandID(mapping), inputLigands(ligands) {}
        double alignmentScore;
        const std::unordered_map<LigandID, PoseID> poseIDsByLigandID;
        std::vector<Ligand> inputLigands;
    };
}  // namespace coaler::multialign
