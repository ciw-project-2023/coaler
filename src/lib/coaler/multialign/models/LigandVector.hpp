#ifndef COALER_LIGANDVECTOR_HPP
#define COALER_LIGANDVECTOR_HPP

#include <vector>

#include "GraphMol/RDKitBase.h"

#include "Ligand.hpp"

namespace coaler::multialign {
    class LigandVector : public std::vector<Ligand> {
      public:
        using std::vector<Ligand>::vector;
        explicit LigandVector(RDKit::MOL_SPTR_VECT molecules);
    };
}
#endif  // COALER_LIGANDVECTOR_HPP
