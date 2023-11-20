//
// Created by chris on 11/20/23.
//

#pragma once
#include <unordered_map>
#include "BasicClasses/Forward.hpp"

namespace MultiAlign
{
    struct MultiAlignerResult{
        double score;
        std::unordered_map<LigandID, PoseID> poseIDsByLigandID;
    };
}