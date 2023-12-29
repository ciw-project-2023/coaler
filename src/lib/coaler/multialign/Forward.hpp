#pragma once
#include <cassert>
#include <unordered_map>

#include "models/Forward.hpp"
#include "models/Ligand.hpp"
#include "models/LigandPair.hpp"
#include "models/PosePair.hpp"

namespace coaler::multialign {

    namespace constants {
        constexpr unsigned DefaultNofStartingAssemblies = 50;
        constexpr unsigned DefaultNofThreads = 1;
        constexpr double PoseRegisterSizeFactor = 0.5;
    }  // namespace constants

    static_assert(constants::PoseRegisterSizeFactor < 1);
    static_assert(constants::PoseRegisterSizeFactor > 0);

    class PoseRegister;
    using PoseRegisterPtr = std::shared_ptr<PoseRegister>;

    using PairwiseAlignment = std::unordered_map<PosePair, double, PosePairHash>;
    using PairwisePoseRegisters = std::unordered_map<LigandPair, PoseRegisterPtr, LigandPairHash>;

    using LigandVector = std::vector<Ligand>;

}  // namespace coaler::multialign
