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
    }  // namespace Constants

    class PoseRegister;
    using PoseRegisterPtr = std::shared_ptr<PoseRegister>;

    using PairwiseAlignment = std::unordered_map<PosePair, double, PosePairHash>;
    using PairwisePoseRegisters = std::unordered_map<LigandPair, PoseRegisterPtr, LigandPairHash>;

    using LigandVector = std::vector<Ligand>;

}  // namespace coaler::multialign
