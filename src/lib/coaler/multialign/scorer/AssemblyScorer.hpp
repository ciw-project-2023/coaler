/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#pragma once
#include "coaler/multialign/Forward.hpp"
#include "coaler/multialign/LigandAlignmentAssembly.hpp"
#include "coaler/multialign/PoseRegisterCollection.hpp"

namespace coaler::multialign {
    class AssemblyScorer {
      public:
        static double calculateAssemblyScore(const coaler::multialign::LigandAlignmentAssembly& assembly,
                                             const coaler::multialign::PairwiseAlignment& scores,
                                             const coaler::multialign::LigandVector& ligands);

        static double calculateScoreDeficitForLigand(LigandID ligandId, LigandID maxLigandId,
                                                     const LigandAlignmentAssembly& assembly,
                                                     const PoseRegisterCollection& registers,
                                                     const PairwiseAlignment& scores);
    };

}  // namespace coaler::multialign
