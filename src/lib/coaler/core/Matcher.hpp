//
// Created by malte on 12/4/23.
//

#pragma once

#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <spdlog/spdlog.h>
#include "GraphMol/RDKitBase.h"
#include "GraphMol/ChemTransforms/ChemTransforms.h"
#include "GraphMol/DistGeomHelpers/Embedder.h"
#include "GraphMol/FMCS/FMCS.h"

namespace coaler::core {
    using AtomMap = std::map<int, RDGeom::Point3D>;
    using CoreMatch = std::vector<RDKit::MatchVectType>;
    using CoreResult = std::pair<RDKit::ROMOL_SPTR, AtomMap>;

    class Matcher {
      public:
        virtual ~Matcher() = default;

        Matcher() = default;

        virtual std::optional<RDKit::MCSResult> findMCS(RDKit::MOL_SPTR_VECT mols);

        std::optional<CoreResult> getResult();

        std::optional<CoreMatch> matchMolecule(RDKit::ROMOL_SPTR mol);

    protected:
        CoreResult core_;
    };
}  // namespace coaler::core