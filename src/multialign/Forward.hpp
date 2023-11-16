//
// Created by chris on 11/4/23.
//

#pragma once
#include <cassert>
#include <boost/functional/hash.hpp>
#include <unordered_map>

namespace MultiAlign
{
    using LigandID = unsigned;
    using PoseID = unsigned;

    /*
    class Ligand;
    class LigandAlignmentAssembly;
    class MultiAligner;
    class PairwiseAlignment;
    class PoseRegisterBuilder; */
    class PosePair;
    struct PosePairHash;

    class PoseRegister;
    using PoseRegisterPtr = std::shared_ptr<PoseRegister>;
    using PairwiseAlignment = std::unordered_map<PosePair, double, PosePairHash>;




}
