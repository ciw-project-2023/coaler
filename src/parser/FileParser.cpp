#include "FileParser.hpp"

#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <spdlog/spdlog.h>

#include <filesystem>
#include <vector>

#include "FileNotFoundException.hpp"

namespace coaler {

    std::vector<RDKit::RWMol *> FileParser::parse(const std::string &file_path) {
        std::vector<RDKit::RWMol *> result;

        spdlog::debug("searching at: {}", std::filesystem::current_path().string());

        std::ifstream infile(file_path);
        if (!infile) {
            spdlog::error("file not found: {}", file_path);

            throw FileNotFoundException(file_path);
        }

        std::string smiles;
        while (std::getline(infile, smiles)) {
            RDKit::RWMol *mol = RDKit::SmilesToMol(smiles);
            RDKit::INT_VECT conf_ids;

            result.push_back(mol);
        }

        return result;
    }
}  // namespace coaler
