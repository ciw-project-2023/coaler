/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */
#pragma once
#include <boost/functional/hash.hpp>
#include <unordered_set>

namespace coaler::multialign {

    class Ligand;

    class LigandPair;

    class PosePair;

    struct UniquePoseID;
    struct UniquePoseIdentifierHash;

    // hash structs
    struct PosePairHash;
    struct LigandPairHash;

    // id declarations
    using LigandID = unsigned;
    using PoseID = unsigned;

    using UniquePoseSet = std::unordered_set<UniquePoseID, UniquePoseIdentifierHash>;

}  // namespace coaler::multialign
