//
// Created by chris on 11/16/23.
//

#pragma once
#include <boost/functional/hash.hpp>
#include <unordered_set>

namespace MultiAlign {

    class Ligand;

    class LigandPair;

    class PosePair;

    struct UniquePoseIdentifier;
    struct UniquePoseIdentifierHash;

    struct PosePairHash;
    struct LigandPairHash;

    using LigandID = unsigned;
    using PoseID = unsigned;

    using UniquePoseSet = std::unordered_set<UniquePoseIdentifier, UniquePoseIdentifierHash>;


}