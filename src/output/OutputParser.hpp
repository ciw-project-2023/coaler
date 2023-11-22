#pragma once

#include <GraphMol/ROMol.h>

namespace coaler {
    /** @struct AlignedMolPair
     *  @brief This structure represents an aligned molecule pair.
     *  @var align_score
     *  @var mol_a
     *  @var mol_b
     */
    struct AlignedMolPair {
        double align_score;
        RDKit::ROMol mol_a;
        RDKit::ROMol mol_b;
    };

    /**
     * @brief This class is responsible for handling the output.
     */
    class OutputParser {
      public:
        /**
         * Add molecule pair with alignment score to aligned molecules
         * @param aligned_mol_pair
         */
        void add_aligned_mols(AlignedMolPair& aligned_mol_pair);

        /**
         * Writes the aligned molecules in an output file.
         * @param file_path
         */
        void save_molecules_w_scores_in_file(const std::string& file_path);

      private:
        std::vector<AlignedMolPair> aligned_mols_;
    };

}  // namespace coaler
