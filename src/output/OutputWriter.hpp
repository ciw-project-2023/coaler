#pragma once

#include <GraphMol/ROMol.h>

namespace coaler {
    /** @struct AlignedMolPair
     *  @brief This structure represents an aligned molecule pair.
     *  @var align_score
     *  @var mol_a
     *  @var mol_b
     *  @var mol_a
     *  @var mol_b
     */
    struct AlignedMolPair {
        double align_score;

        int id_mol_a;
        int id_mol_b;

        RDKit::ROMol mol_a;
        RDKit::ROMol mol_b;
    };

    /**
     * @brief This class is responsible for handling the output.
     */
    class OutputWriter {
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

        /**
         * Prints the aligned molecules with score inside the log.
         */
        void print_in_log_molecules_w_scores();

      private:
        std::vector<AlignedMolPair> aligned_mols_;
    };

}  // namespace coaler
