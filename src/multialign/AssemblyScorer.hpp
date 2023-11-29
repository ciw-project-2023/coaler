//
// Created by chris on 11/27/23.
//

#pragma once
#include "Forward.hpp"
#include "LigandAlignmentAssembly.hpp"
#include "PoseRegisterCollection.hpp"

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
