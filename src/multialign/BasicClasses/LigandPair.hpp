//
// Created by chris on 11/9/23.
//

#pragma once
#include "Forward.hpp"

    namespace coaler::multialign
    {
        /**
         * A pair of ligands.
         */
        class LigandPair{
        public:
            LigandPair(LigandID first,
                       LigandID second);

            [[nodiscard]] LigandID getFirst() const noexcept;
            [[nodiscard]] LigandID getSecond() const noexcept;

            bool operator==(const LigandPair& other) const;

        private:
            LigandID m_firstLigand;
            LigandID m_secondLigand;
        };

        struct LigandPairHash
        {
            std::size_t operator()(const LigandPair& pair) const
            {
                std::size_t seed = 0;

                boost::hash_combine(seed,boost::hash_value(pair.getFirst()));
                boost::hash_combine(seed,boost::hash_value(pair.getSecond()));

                return seed;
            }};
}
