#include "OutputParser.hpp"

#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <spdlog/spdlog.h>

namespace coaler {
    void OutputParser::add_aligned_mols(AlignedMolPair& aligned_mol_pair) {
        aligned_mols_.emplace_back(aligned_mol_pair);
    }

    void OutputParser::save_molecules_w_scores_in_file(const std::string& file_path) {
        // TODO: define here output format
    }

    void OutputParser::print_in_log_molecules_w_scores() {
        spdlog::info("Start printing results.");
        for (auto mol_pair : aligned_mols_) {
            spdlog::info("Molecule {}:{} and molecule {}:{} are aligned with score: {}", mol_pair.id_mol_a,
                         MolToSmarts(mol_pair.mol_a), mol_pair.id_mol_b, MolToSmarts(mol_pair.mol_b),
                         mol_pair.align_score);
        }
        spdlog::info("End of results. :)");
    }
}  // namespace coaler
