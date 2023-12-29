#pragma once

#include <GraphMol/ROMol.h>
#include <GraphMol/Substruct/SubstructMatch.h>

namespace coaler::multialign {

    /**
     * @brief This class is responsible the scoring of the alignment of two molecules.
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
        static double calc_tanimoto_shape_similarity(const RDKit::ROMol& mol_a, const RDKit::ROMol& mol_b,
                                                     unsigned int pos_id_a, unsigned int pos_id_b);

        static double calc_alignment_transformation(const RDKit::ROMol& mol_a, const RDKit::ROMol& mol_b,
                                                    unsigned int pos_id_a, unsigned int pos_id_b);
        static double calc_rmsd(RDKit::ROMol mol_a, RDKit::ROMol mol_b, int pos_id_a, int pos_id_b);
    };

}  // namespace coaler::multialign