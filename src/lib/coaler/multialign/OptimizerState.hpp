#pragma once

#include "LigandAlignmentAssembly.hpp"
#include "PoseRegisterCollection.hpp"
#include "models/Forward.hpp"

namespace coaler::multialign {

    /**
     * The OptimizerState struct represents an assembly along with everything required for the optimization.
     */
    // NOLINTNEXTLINE(cppcoreguidelines-special-member-functions,)
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