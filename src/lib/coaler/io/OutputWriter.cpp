/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#include "OutputWriter.hpp"

#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/RWMol.h>

#include <spdlog/spdlog.h>


namespace coaler::io {

    void OutputWriter::save_molecules_w_scores_in_file(const std::string &file_path, const coaler::multialign::MultiAlignerResult& result) {
        std::ofstream output_file(file_path);
        if (!output_file.is_open()) {
            spdlog::error("Cannot open file: {}", file_path);
            return;
        }

        boost::shared_ptr<RDKit::SDWriter> const sdf_writer(new RDKit::SDWriter(&output_file, false));
        if (result.poseIDsByLigandID.size() != result.inputLigands.size()) {
            spdlog::info("only generated an incomplete alignment.");
            return;
        }

        for (const auto &entry: result.inputLigands) {
            sdf_writer->write(entry.getMolecule(), result.poseIDsByLigandID.at(entry.getID()));
        }
    }
}  // namespace coaler::io
