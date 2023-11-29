//
// Created by chris on 11/27/23.
//

#pragma once
#include "Forward.hpp"
#include "LigandAlignmentAssembly.hpp"
#include "PoseRegisterCollection.hpp"


namespace MultiAlign
{
class AssemblyScorer {

    static double calculateAssemblyScore(
        const MultiAlign::LigandAlignmentAssembly& assembly,
        const MultiAlign::PairwiseAlignment& scores,
        const MultiAlign::LigandVector& ligands
    );

    static double calculateScoreDeficitForLigand(
        LigandID ligandId,
        LigandID maxLigandId,
        const LigandAlignmentAssembly& assembly,
        const PoseRegisterCollection& registers,
        const PairwiseAlignment& scores
        );

};

}
