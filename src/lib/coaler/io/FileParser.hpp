/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <spdlog/spdlog.h>

#include <filesystem>

#include "FileNotFoundException.hpp"

namespace coaler::io {
    class FileParser {
      public:
        static std::vector<RDKit::RWMol *> parse(const std::string &file_path);
    };
}  // namespace coaler::io
