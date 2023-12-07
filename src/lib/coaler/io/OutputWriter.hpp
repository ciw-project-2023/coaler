/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#pragma once

#include <GraphMol/ROMol.h>

#include <unordered_map>

#include "../multialign/Forward.hpp"
#include "../multialign/MultiAlignerResult.hpp"

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
        static void writeSDF(const std::string &file_path, const coaler::multialign::MultiAlignerResult &result);
    };

}  // namespace coaler::io
