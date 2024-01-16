#pragma once

#include "models/Forward.hpp"
#include "LigandAlignmentAssembly.hpp"
#include "PoseRegisterCollection.hpp"

namespace coaler::multialign {

    /**
     * The OptimizerState struct represents an assembly along with everything required for the optimization.
     */
    struct OptimizerState {
        // NOLINTBEGIN(misc-non-private-member-variables-in-classes)
        double score;
        LigandAlignmentAssembly assembly;
        PairwiseAlignments scores;
        LigandVector ligands;
        PoseRegisterCollection registers;
        // NOLINTEND(misc-non-private-member-variables-in-classes)

        OptimizerState& operator=(const OptimizerState& s) = default;
    };
}  // namespace coaler::multialign