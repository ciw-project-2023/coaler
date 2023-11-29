/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#include "OutputWriter.hpp"

#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/RWMol.h>
#include <spdlog/spdlog.h>

namespace coaler::io {
    void OutputWriter::save_molecules_w_scores_in_file(MultiAlignerResult result, const std::string& file_path) {
#pragma unroll 100
        for (auto ligand_id_conformer_pair : result.poseIDsByLigandID) {
            LigandID ligand_id = std::get<0>(ligand_id_conformer_pair);
            PoseID poseId = std::get<1>(ligand_id_conformer_pair);
            RDKit::RWMol* ligand_w_ptr = result.ligands.at(ligand_id);
            const RDKit::ROMol ligand_o = *ligand_w_ptr;
        }

        // TODO: write sdf file
    }

    void OutputWriter::print_multi_aligner_result(MultiAlignerResult result) {
        spdlog::info("Start printing results.");

#pragma unroll 100
        for (auto ligand_id_conformer_pair : result.poseIDsByLigandID) {
            LigandID ligand_id = std::get<0>(ligand_id_conformer_pair);
            PoseID poseId = std::get<1>(ligand_id_conformer_pair);
            RDKit::RWMol* ligand_w_ptr = result.ligands.at(ligand_id);
            const RDKit::ROMol ligand_o = *ligand_w_ptr;

            spdlog::info("Molecule Smiles:{} ID:{} Conformere: {}", RDKit::MolToSmiles(ligand_o), ligand_id, poseId);
        }
        spdlog::info("Multi Alginement Score: {}", result.score);
        spdlog::info("End of results. :)");
    }
}  // namespace coaler::io
