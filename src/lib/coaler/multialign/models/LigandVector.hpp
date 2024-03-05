#pragma once
#include <vector>

#include "GraphMol/RDKitBase.h"
#include "Ligand.hpp"

namespace coaler::multialign {

     /**
      * @brief LigandVector class to represent a vector of ligands.
      *
      * The LigandVector class is used to represent a vector of ligands.
      * It is used to store the ligands in a single object.
      */
    class LigandVector : public std::vector<Ligand> {
      public:
        using std::vector<Ligand>::vector;
        explicit LigandVector(RDKit::MOL_SPTR_VECT& molecules);
    };
}  // namespace coaler::multialign
