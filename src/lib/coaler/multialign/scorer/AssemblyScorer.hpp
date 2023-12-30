/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#pragma once
#include <coaler/multialign/models/PairwiseAlignments.hpp>

#include "coaler/multialign/Forward.hpp"
#include "coaler/multialign/LigandAlignmentAssembly.hpp"
#include "coaler/multialign/PoseRegisterCollection.hpp"

namespace coaler::multialign {
    class AssemblyScorer {
      public:
        static double calculateAssemblyScore(const coaler::multialign::LigandAlignmentAssembly& assembly,
                                             coaler::multialign::PairwiseAlignments& scores,
                                             const coaler::multialign::LigandVector& ligands);

        static double calculateScoreDeficitForLigand(LigandID ligandId, LigandID maxLigandId,
                                                     const LigandAlignmentAssembly& assembly,
                                                     const PoseRegisterCollection& registers,
                                                     PairwiseAlignments& scores, const LigandVector& ligands);

        static double getScoreInAssembly(LigandID firstLigandID, LigandID secondLigandID, PoseID firstPoseID,
                                         PoseID secondPoseID, PairwiseAlignments& scores, const LigandVector& ligands);
    };

}  // namespace coaler::multialign
