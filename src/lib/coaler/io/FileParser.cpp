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

    RDKit::MOL_SPTR_VECT FileParser::parse(const std::string& filePath) {
        RDKit::MOL_SPTR_VECT result;

        spdlog::debug("searching at: {}", std::filesystem::current_path().string());

        const std::ifstream infile(filePath);
        if (!infile) {
            spdlog::error("file not found: {}", filePath);
            throw FileNotFoundException(filePath);
        }

        const auto fileExtension = std::filesystem::path(filePath).extension().string();
        if (fileExtension == ".sdf") {
            // Parse as sdf file
            RDKit::SDMolSupplier supplier(filePath, true, false, false);
            while (!supplier.atEnd()) {
                auto mol = boost::make_shared<RDKit::RWMol>(*supplier.next());
                result.emplace_back(mol);
            }
        } else if (fileExtension == ".smi") {
            // Parse as smi file
            int line = 0;
            RDKit::SmilesMolSupplier supplier(filePath, "\t", 0, 1, false, true);
            while (!supplier.atEnd()) {
                line++;

                auto* const next = supplier.next();
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

        return checkInputMolecules(result, filePath);
    }

    RDKit::MOL_SPTR_VECT FileParser::checkInputMolecules(const RDKit::MOL_SPTR_VECT& mols,
                                                         const std::string& filePath) {
        RDKit::MOL_SPTR_VECT retMols;
        std::unordered_set<std::string> smilesSet;
        unsigned duplicates = 0;
        for (auto const& mol : mols) {
            std::string smiles = RDKit::MolToSmiles(*mol);
            if (smilesSet.insert(smiles).second) {
                retMols.emplace_back(mol);
            } else {
                duplicates++;
            }
        }

        if (duplicates > 0) {
            spdlog::warn("input file {} contains {} duplicates. All duplicates will be removed", filePath, duplicates);
        }

        return retMols;
    }
}  // namespace coaler::io
