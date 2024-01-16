#pragma once
#include <vector>

#include "GraphMol/RDKitBase.h"
#include "Ligand.hpp"

namespace coaler::multialign {
    class LigandVector : public std::vector<Ligand> {
      public:
        using std::vector<Ligand>::vector;
        explicit LigandVector(RDKit::MOL_SPTR_VECT& molecules);
    };
}  // namespace coaler::multialign
