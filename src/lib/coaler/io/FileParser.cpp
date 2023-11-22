#include "FileParser.hpp"

#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <spdlog/spdlog.h>

#include <filesystem>
#include <vector>

#include "FileNotFoundException.hpp"

namespace coaler::io {

    std::vector<RDKit::RWMol *> FileParser::parse(const std::string &file_path) {
        std::vector<RDKit::RWMol *> result;

        spdlog::debug("searching at: {}", std::filesystem::current_path().string());

        std::ifstream infile(file_path);
        if (!infile) {
            spdlog::error("file not found: {}", file_path);

            throw FileNotFoundException(file_path);
        }

        std::string file_extension = std::filesystem::path(file_path).extension().string();
        if (file_extension == ".sdf") {
            // Parse as sdf file
            RDKit::SDMolSupplier supplier(file_path, true, false, false);
            while (!supplier.atEnd()) {
                RDKit::RWMol *mol = new RDKit::RWMol(*supplier.next());
                result.push_back(mol);
            }
        } else if (file_extension == ".smi") {
            // Parse as smi file
            RDKit::SmilesMolSupplier supplier(file_path, "\t", 0, -1, false, true);
            while (!supplier.atEnd()) {
                RDKit::RWMol *mol = new RDKit::RWMol(*supplier.next());
                result.push_back(mol);
            }
        } else {
            throw std::runtime_error("Unsupported file extension: " + file_extension);
        }

        return result;
    }
}  // namespace coaler::io
