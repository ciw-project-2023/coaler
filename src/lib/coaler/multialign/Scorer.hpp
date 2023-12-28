//
// Created by cnoack on 28.12.23.
//
#include <GraphMol/ShapeHelpers/ShapeUtils.h>

#include "Forward.hpp"
#include "coaler/multialign/models/Ligand.hpp"

namespace coaler::multialign {
    class Scorer {
      public:
        static inline double getOverlapScore(const Ligand &ligand1, const Ligand &ligand2, unsigned int pose1,
                                             unsigned int pose2) {
            const double distance
                = RDKit::MolShapes::tanimotoDistance(ligand1.getMolecule(), ligand2.getMolecule(), pose1, pose2);
            const double similarity = 1 - distance;
            return similarity;
        }
    };

}  // namespace coaler::multialign
