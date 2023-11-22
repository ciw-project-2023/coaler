#pragma once

#include <GraphMol/ROMol.h>

namespace coaler {
    struct AlignedMolPair {
        double align_score;
        RDKit::ROMol mol_a;
        RDKit::ROMol mol_b;
    };

    class OutputParser {
      public:
        void add_aligned_mols(AlignedMolPair& aligned_mol_pair);

        void save_molecules_w_scores_in_file(const std::string& file_path);

      private:
        std::vector<AlignedMolPair> aligned_mols_;
    };

}  // namespace coaler
