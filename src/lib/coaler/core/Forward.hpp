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

namespace coaler::core {

    class Core;

    enum class coreType {
        MCS [[maybe_unused]],
        Murcko [[maybe_unused]],
        MurckoGeneric [[maybe_unused]],
    };

    using CoreAsMol = RDKit::ROMol;
    using CorePtr = RDKit::ROMOL_SPTR;
    using CoreMolMatch = std::vector<RDKit::MatchVectType>;

}  // namespace coaler::core