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
        static RDKit::MOL_SPTR_VECT parse(const std::string &file_path);
      private:
        /***
         * checks the input molecules for duplicates and deletes them
         * @param mols molecules parsed in parse()
         * @return RDKit::MOL_SPTR_VECT without duplicates
         */
        static RDKit::MOL_SPTR_VECT checkInputMolecules(const RDKit::MOL_SPTR_VECT mols);
    };
}  // namespace coaler::io
