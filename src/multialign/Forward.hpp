//
// Created by chris on 11/4/23.
//

#pragma once
#include <cassert>
#include <boost/functional/hash.hpp>
#include "boost/numeric/ublas/triangular.hpp"
//#include "LigandAlignmentAssembly.hpp"
#include "PairwiseAlignment.hpp"

namespace MultiAlign
{
    using LigandID = unsigned;
    using PoseID = unsigned;

    /*
    class Ligand;
    class LigandAlignmentAssembly;
    class MultiAligner;
    class PairwiseAlignment;
    class PosePair;
    class PoseRegisterBuilder; */

    class PoseRegister;
    using PoseRegisterPtr = std::shared_ptr<PoseRegister>;




}
