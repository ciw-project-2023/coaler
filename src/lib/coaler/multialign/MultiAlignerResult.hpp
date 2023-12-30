/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

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

            /*--------------------------------------------------------------------------------------------------------*/

            MultiAlignerResult& operator=(const MultiAlignerResult& result) {
                this->alignmentScore = result.alignmentScore;
                this->poseIDsByLigandID = result.poseIDsByLigandID;
                this->inputLigands = result.inputLigands;
                return *this;
            };

            /*--------------------------------------------------------------------------------------------------------*/

            double alignmentScore;
            std::unordered_map<LigandID, PoseID> poseIDsByLigandID;
            std::vector<Ligand> inputLigands;
        };

    }  // namespace multialign
}  // namespace coaler
