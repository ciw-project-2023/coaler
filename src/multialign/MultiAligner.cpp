//
// Created by chris on 11/5/23.
//

#include "MultiAligner.hpp"

namespace MultiAlign
{
    MultiAligner::MultiAligner(
            const std::vector<RDKit::RWMol>& molecules,
            const RDKit::MCSResult& core)
            : m_molecules(molecules)
            , m_core(core)
            {
            }

    /*----------------------------------------------------------------------------------------------------------------*/

    MultiAlignerResult MultiAligner::alignMolecules(){
        //PairwiseAlignment allPosesAlignmentScores = SingleAligner(...);

        //generate registers

        // build starting ensembles from registers

        //optimize ensembles

        return MultiAlignerResult();
    }
}