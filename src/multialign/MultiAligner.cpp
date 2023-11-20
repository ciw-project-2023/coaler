//
// Created by chris on 11/5/23.
//

#include "MultiAligner.hpp"

namespace MultiAlign
{
    MultiAligner::MultiAligner(
            const std::vector<RDKit::RWMol>& molecules,
            const RDKit::MCSResult& core) {

    }

    /*----------------------------------------------------------------------------------------------------------------*/

    MultiAlignerResult MultiAligner::alignMolecules() {
        return MultiAlignerResult();
    }
}