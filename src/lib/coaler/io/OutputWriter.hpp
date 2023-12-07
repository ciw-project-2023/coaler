/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#pragma once

#include <GraphMol/ROMol.h>
#include "../multialign/Forward.hpp"
#include "../multialign/MultiAlignerResult.hpp"

#include <unordered_map>

namespace coaler::io {
    /**
     * @brief This class is responsible for handling the output.
     */
    class OutputWriter {
      public:
        /**
         * Writes the aligned molecules in an output file.
         * @param file_path
         */
        static void save_molecules_w_scores_in_file(const std::string &file_path, const multialign::MultiAlignerResult &result);
    };

}  // namespace coaler::io
