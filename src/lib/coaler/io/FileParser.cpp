/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#include "FileParser.hpp"

#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <spdlog/spdlog.h>

#include <filesystem>
#include <vector>

#include "FileNotFoundException.hpp"

namespace coaler::io {

    RDKit::MOL_SPTR_VECT FileParser::parse(const std::string &file_path) {
        RDKit::MOL_SPTR_VECT result;

        spdlog::debug("searching at: {}", std::filesystem::current_path().string());

        std::ifstream infile(file_path);
        if (!infile.is_open()) {
            spdlog::error("file not found: {}", file_path);
            throw FileNotFoundException(file_path);
        }

        const std::string fileExtension = std::filesystem::path(file_path).extension().string();
        if (fileExtension == ".sdf") {
            // Parse as sdf file
            RDKit::SDMolSupplier supplier(file_path, true, false, false);
            while (!supplier.atEnd()) {
                auto mol = boost::make_shared<RDKit::RWMol>(*supplier.next());
                result.emplace_back(mol);
            }
        } else if (fileExtension == ".smi") {
            // Parse as smi file
            int line = 0;
            RDKit::SmilesMolSupplier supplier(file_path, "\t", 0, 1, false, true);
            while (!supplier.atEnd()) {
                line++;

                auto *const next = supplier.next();
                if (next == nullptr) {
                    spdlog::warn("Could not parse line {} in smi file", line);
                    continue;
                }

                auto mol = boost::make_shared<RDKit::RWMol>(*next);
                result.emplace_back(mol);
            }
        } else {
            spdlog::error("Unsupported file extension: {}", fileExtension);
            throw std::runtime_error("Unsupported file extension: " + fileExtension);
        }

        return result;
    }
}  // namespace coaler::io
