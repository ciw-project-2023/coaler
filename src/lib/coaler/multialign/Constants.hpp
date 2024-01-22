#pragma once

namespace coaler::multialign::constants {
    const unsigned DEFAULT_NOF_STARTING_ASSEMBLIES = 50;
    const unsigned DEFAULT_NOF_THREADS = 1;
    const double POSE_REGISTER_SIZE_FACTOR = 0.5;

    /**
     * This treshold determines the score deficit above which new conformers are attempted to be generated
     * during assembly optimization.
     */
    const double COARSE_OPTIMIZATION_THRESHOLD = 1.5;

    /**
     * This treshold determines the score deficit above which new conformers are attempted to be generated
     * during the fine tuning of the best alignment assembly.
     */
    const double FINE_OPTIMIZATION_THRESHOLD = 0.1;

    const unsigned OPTIMIZER_STEP_LIMIT = 100;
}  // namespace coaler::multialign::constants
