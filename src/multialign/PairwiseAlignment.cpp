//
// Created by chris on 11/4/23.
//

#include "PairwiseAlignment.hpp"


namespace MultiAlign
{

    PairwiseAlignment::PairwiseAlignment(const PairwiseAlignmentMatrix& singleAlignmentScores)
    : m_singleAlignmentScores(singleAlignmentScores)
    {}

    /*---------------------------------------------------------------------------------*/

    double PairwiseAlignment::getValue(const unsigned int x, const unsigned int y) const {
        assert(x < m_singleAlignmentScores.size1()
               && y < m_singleAlignmentScores.size2());
        assert(y > x);
        return m_singleAlignmentScores(x,y);
    }
}