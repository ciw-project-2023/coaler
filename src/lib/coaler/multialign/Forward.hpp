/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#pragma once
#include <cassert>
#include <unordered_map>

#include "models/Forward.hpp"
#include "models/Ligand.hpp"
#include "models/LigandPair.hpp"
#include "models/PosePair.hpp"

namespace coaler::multialign {

    namespace Constants {
        constexpr unsigned DEFAULT_NOF_STARTING_ASSEMBLIES = 50;
        constexpr unsigned DEFAULT_NOF_THREADS = 1;
        constexpr double POSE_REGISTER_SIZE_FACTOR = 0.5;
    }  // namespace Constants

    static_assert(Constants::POSE_REGISTER_SIZE_FACTOR < 1);
    static_assert(Constants::POSE_REGISTER_SIZE_FACTOR > 0);

    class PoseRegister;
    using PoseRegisterPtr = std::shared_ptr<PoseRegister>;

    using PairwiseAlignment = std::unordered_map<PosePair, double, PosePairHash>;
    using PairwisePoseRegisters = std::unordered_map<LigandPair, PoseRegisterPtr, LigandPairHash>;

    using LigandPtr = boost::shared_ptr<Ligand>;
    using LigandVector = std::vector<Ligand>;

}  // namespace coaler::multialign
