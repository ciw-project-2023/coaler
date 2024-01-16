#pragma once

#include <unordered_map>
#include <vector>

#include "Alias.hpp"
#include "LigandVector.hpp"
#include "PosePair.hpp"

namespace coaler::multialign {

    /**
     * This class stores pairwise conformer overlap values. If its presented a pair that it hasn´t encounted
     * before, it calculates the overlap and stores it if desired.
     */
    // NOLINTBEGIN(cppcoreguidelines-virtual-class-destructor, cppcoreguidelines-special-member-functions)
    class PairwiseAlignments : public std::unordered_map<PosePair, double, PosePairHash> {
      public:
        PairwiseAlignments() = default;
        ~PairwiseAlignments() = default;
        PairwiseAlignments(PairwiseAlignments& p);
        PairwiseAlignments(const PairwiseAlignments& p);
        PairwiseAlignments(PairwiseAlignments&& p) = default;

        /**
         * @brief looks up or calculates the overlap score
         *
         * @param key The pair to score
         * @param ligands The ligands
         * @param store
         * @return
         */
        double at(const PosePair& key, const LigandVector& ligands = {}, bool store = false);

        PairwiseAlignments& operator=(const PairwiseAlignments& p);
        // NOLINTNEXTLINE(readability-avoid-const-params-in-decls)
        virtual PairwiseAlignments& operator=(const std::unordered_map<PosePair, double, PosePairHash>& p);
    };
    // NOLINTEND(cppcoreguidelines-virtual-class-destructor, cppcoreguidelines-special-member-functions)
}  // namespace coaler::multialign
