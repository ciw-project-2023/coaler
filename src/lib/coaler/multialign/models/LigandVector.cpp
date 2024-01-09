//
// Created by niklas on 1/9/24.
//

#include "LigandVector.hpp"
#include "Alias.hpp"
#include "UniquePoseSet.hpp"

namespace  coaler::multialign {
    LigandVector::LigandVector(RDKit::MOL_SPTR_VECT molecules) {
        for (LigandID id = 0; id < molecules.size(); id++) {
            UniquePoseSet poses;

            for (PoseID poseId = 0; poseId < molecules.at(id)->getNumConformers(); poseId++) {
                poses.emplace(id, poseId);
            }

            auto mol = *molecules.at(id);
            const Ligand ligand(mol, poses, id);

            this->emplace_back(ligand);
        }
    }
}
