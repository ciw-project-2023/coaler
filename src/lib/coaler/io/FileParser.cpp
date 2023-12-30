/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#include "FileParser.hpp"

#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <spdlog/spdlog.h>

#include <filesystem>
#include <vector>

#include "FileNotFoundException.hpp"

namespace coaler::io {

    RDKit::MOL_SPTR_VECT FileParser::parse(const std::string &file_path) {
        RDKit::MOL_SPTR_VECT result;

        spdlog::debug("searching at: {}", std::filesystem::current_path().string());

        std::ifstream infile(file_path);
        if (!infile) {
            spdlog::error("file not found: {}", file_path);
            throw FileNotFoundException(file_path);
        }

        const auto file_extension = std::filesystem::path(file_path).extension().string();
        if (file_extension == ".sdf") {
            // Parse as sdf file
            RDKit::SDMolSupplier supplier(file_path, true, false, false);
            while (!supplier.atEnd()) {
                auto mol = boost::make_shared<RDKit::RWMol>(*supplier.next());
                result.emplace_back(mol);
            }
        } else if (file_extension == ".smi") {
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
            spdlog::error("Unsupported file extension: {}", file_extension);
            throw std::runtime_error("Unsupported file extension: " + file_extension);
        }

        return checkInputMolecules(result);
    }
    RDKit::MOL_SPTR_VECT FileParser::checkInputMolecules(const RDKit::MOL_SPTR_VECT mols) {

        RDKit::MOL_SPTR_VECT retMols;
        std::vector<std::string> smilesVec;
        unsigned duplicates;
        for (auto mol : mols) {
            std::string smiles = RDKit::MolToSmiles(*mol);
            if (std::find(smilesVec.begin(), smilesVec.end(), smiles) == smilesVec.end()) {
                smilesVec.emplace_back(smiles);
                retMols.emplace_back(mol);
            }
            else {
                duplicates++;
            }
        }
        if (duplicates > 0) {
            spdlog::warn("input file contains {} duplicates. All duplicates will be removed", duplicates);
        }

        return retMols;

    }
}  // namespace coaler::io
