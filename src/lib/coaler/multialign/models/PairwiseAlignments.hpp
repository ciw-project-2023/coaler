//
// Created by chris on 12/29/23.
//

#pragma once

#include <unordered_map>
#include <vector>

#include "Forward.hpp"
#include "PosePair.hpp"
namespace coaler::multialign{

class PairwiseAlignments : public std::unordered_map<PosePair, double, PosePairHash>{

  public:

    PairwiseAlignments() = default;
    PairwiseAlignments(PairwiseAlignments& p);
    PairwiseAlignments(const PairwiseAlignments& p);


    double at(const PosePair& key, const std::vector<Ligand>& ligands = {}, bool store = false);

    PairwiseAlignments& operator= (const PairwiseAlignments& p);
    virtual PairwiseAlignments& operator= (const std::unordered_map<PosePair, double, PosePairHash>& p);
};
}

