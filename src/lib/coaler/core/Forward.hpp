//
// Created by malte on 12/4/23.
//

#pragma once

#include <cassert>
#include <string>
#include <vector>
#include <utility>
#include <GraphMol/RWMol.h>

namespace coaler::core {

    class Core;

    enum class coreType {
        MCS [[maybe_unused]],
        Murcko [[maybe_unused]],
        MurckoGeneric [[maybe_unused]],
    };

    using CoreAsMol = RDKit::RWMol;

}