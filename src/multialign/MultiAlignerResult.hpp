//
// Created by chris on 11/20/23.
//

#pragma once
#include <unordered_map>
#include "BasicClasses/Forward.hpp"

namespace coaler{
namespace multialign{

    struct MultiAlignerResult{
        double score;
        std::unordered_map<LigandID, PoseID> poseIDsByLigandID;
        std::vector<Ligand> ligands;
    };

}//multialign
}//coaler