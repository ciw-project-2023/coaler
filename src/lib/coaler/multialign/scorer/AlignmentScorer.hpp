#pragma once

#include <GraphMol/ROMol.h>
#include <GraphMol/Substruct/SubstructMatch.h>

namespace coaler::multialign {

    /**
     * @brief AlignmentScorer class to calculate the tanimoto shape similarity between two molecules
     */
    class AlignmentScorer {
      public:
        /**
         * Computes the tanimoto shape similarity between two molecules.
         * @param mol_a
         * @param mol_b
         * @param pos_id_a: Conforemere ID for molecules A
         * @param pos_id_b: Conforemere ID for molecules B
         * @return Tanimoto shape similarity of molecules
         */
        static double calcTanimotoShapeSimilarity(const RDKit::ROMol& molA, const RDKit::ROMol& molB,
                                                  unsigned int posIdA, unsigned int posIdB);
    };

}  // namespace coaler::multialign