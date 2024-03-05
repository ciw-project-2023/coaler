#pragma once
#include <cstddef>

#include "Alias.hpp"
#include "boost/functional/hash.hpp"

namespace coaler::multialign {
    /*!
     * @file LigandPair.hpp
     * @brief This file contains the LigandPair class which is used to represent a pair of ligands.
     */
    class LigandPair {
      public:
        LigandPair(LigandID first, LigandID second);

        [[nodiscard]] LigandID getFirst() const noexcept;
        [[nodiscard]] LigandID getSecond() const noexcept;

        bool operator==(const LigandPair& other) const;

      private:
        LigandID m_firstLigand;
        LigandID m_secondLigand;
    };
    /**
     * @brief Hash function for LigandPair
     */
    struct LigandPairHash {
        std::size_t operator()(const LigandPair& pair) const {
            std::size_t seed = 0;

            boost::hash_combine(seed, boost::hash_value(pair.getFirst()));
            boost::hash_combine(seed, boost::hash_value(pair.getSecond()));

            return seed;
        }
    };
}  // namespace coaler::multialign
