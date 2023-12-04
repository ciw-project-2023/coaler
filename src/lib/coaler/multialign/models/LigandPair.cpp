/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#include "LigandPair.hpp"

namespace coaler::multialign {
    LigandPair::LigandPair(const LigandID first, const LigandID second) {
        // swap values to rule out duplicates
        assert(first != second);
        if (first > second) {
            m_firstLigand = second;
            m_secondLigand = first;
        } else {
            m_firstLigand = first;
            m_secondLigand = second;
        }
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    bool LigandPair::operator==(const LigandPair &other) const {
        return this->m_firstLigand == other.m_firstLigand && this->m_secondLigand == other.m_secondLigand;
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    LigandID LigandPair::getFirst() const noexcept { return m_firstLigand; }

    /*----------------------------------------------------------------------------------------------------------------*/

    LigandID LigandPair::getSecond() const noexcept { return m_secondLigand; }
}  // namespace coaler::multialign
