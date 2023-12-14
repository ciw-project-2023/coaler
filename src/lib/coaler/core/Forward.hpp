//
// Created by malte on 12/4/23.
//

#pragma once

#include <GraphMol/RWMol.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <cassert>
#include <string>
#include <utility>
#include <vector>

#include "Matcher.hpp"

namespace coaler::core {
    enum class CoreType {
        MCS [[maybe_unused]],
        Murcko [[maybe_unused]],
        MurckoGeneric [[maybe_unused]],
    };

    using CorePtr = RDKit::ROMOL_SPTR;
    using CoreAsMol = CorePtr;
    using CoreMolMatch = std::vector<RDKit::MatchVectType>;

}  // namespace coaler::core