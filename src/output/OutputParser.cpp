#include "OutputParser.hpp"

namespace coaler {
    void OutputParser::add_aligned_mols(AlignedMolPair& aligned_mol_pair) {
        aligned_mols_.emplace_back(aligned_mol_pair);
    }

    void OutputParser::save_molecules_w_scores_in_file(const std::string& file_path) {
        // TODO: define here output format
    }
}  // namespace coaler
