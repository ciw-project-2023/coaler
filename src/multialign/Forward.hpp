//
// Created by chris on 11/4/23.
//

#pragma once
#include <cassert>
#include <unordered_map>

#include "BasicClasses/Forward.hpp"
#include "BasicClasses/Ligand.hpp"
#include "BasicClasses/LigandPair.hpp"
#include "BasicClasses/PosePair.hpp"

namespace MultiAlign
{


    class PoseRegister;
    using PoseRegisterPtr = std::shared_ptr<PoseRegister>;

    using PairwiseAlignment = std::unordered_map<PosePair, double, PosePairHash>;
    using PairwisePoseRegisters = std::unordered_map<LigandPair, PoseRegisterPtr, LigandPairHash>;





}
