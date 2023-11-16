//
// Created by chris on 11/16/23.
//

#pragma once
#include <boost/functional/hash.hpp>

namespace MultiAlign {

    class Ligand;

    class LigandPair;

    class PosePair;

    struct PosePairHash;
    struct LigandPairHash;


    using PoseID = unsigned;
    using LigandID = unsigned;

}