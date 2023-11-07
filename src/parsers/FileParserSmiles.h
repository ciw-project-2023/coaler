//
// Created by niklas on 11/7/23.
//

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <spdlog/spdlog.h>
#include <filesystem>
#include "FileNotFoundException.h"

class FileParserSmiles {
public:
    static std::vector<RDKit::RWMol*> parse(const std::string& file_path);
};

