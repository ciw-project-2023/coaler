/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#pragma once

#include <GraphMol/ROMol.h>

#include <unordered_map>

namespace coaler::io {

    // TODO: Remove after MultiAlign is merged
    using LigandID = unsigned;
    using PoseID = unsigned;

    struct MultiAlignerResult{
        double score;

        std::unordered_map<LigandID, PoseID> poseIDsByLigandID;
        std::vector<RDKit::RWMol *> ligands;
    };
    // TODO: END

    /**
     * @brief This class is responsible for handling the output.
     */
    class OutputWriter {
      public:
        /**
         * Writes the aligned molecules in an output file.
         * @param file_path
         */
        static void save_molecules_w_scores_in_file(MultiAlignerResult result, const std::string& file_path);

        /**
         * Prints the aligned molecules with score inside the log.
         */
        static void print_multi_aligner_result(MultiAlignerResult result);
    };

}  // namespace coaler::io
