//
// Created by chris on 12/29/23.
//

#pragma once

#include <unordered_map>
#include <vector>

#include "Forward.hpp"
#include "PosePair.hpp"
namespace coaler::multialign {

    /**
     * This class stores pairwise conformer overlap values. If its presented a pair that it hasn´t encounted
     * before, it calculates the overlap and stores it if desired.
     */
    class PairwiseAlignments : public std::unordered_map<PosePair, double, PosePairHash> {
      public:
        PairwiseAlignments() = default;
        ~PairwiseAlignments() = default;
        PairwiseAlignments(PairwiseAlignments& p);
        PairwiseAlignments(const PairwiseAlignments& p);

        /**
         * @brief looks up or calculates the overlap score
         *
         * @param key The pair to score
         * @param ligands The ligands
         * @param store
         * @return
         */
        double at(const PosePair& key, const std::vector<Ligand>& ligands = {}, bool store = false);

        PairwiseAlignments& operator=(const PairwiseAlignments& p);
        virtual PairwiseAlignments& operator=(const std::unordered_map<PosePair, double, PosePairHash>& p);
    };
}  // namespace coaler::multialign
