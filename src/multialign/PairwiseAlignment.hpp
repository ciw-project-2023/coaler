//
// Created by chris on 11/4/23.
//

#pragma once
#include <boost/numeric/ublas/triangular.hpp>


namespace MultiAlign
{
    using PairwiseAlignmentMatrix =
    boost::numeric::ublas::triangular_matrix<double,
            boost::numeric::ublas::unit_upper,
            boost::numeric::ublas::column_major>;

    class PairwiseAlignment {

    public:
        PairwiseAlignment(const PairwiseAlignmentMatrix & singleAlignmentScores);

        double getValue(const unsigned x, const unsigned y) const;


    private:
        PairwiseAlignmentMatrix m_singleAlignmentScores;
    };

}


