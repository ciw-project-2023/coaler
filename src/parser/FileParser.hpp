#include <GraphMol/SmilesParse/SmilesParse.h>
#include <spdlog/spdlog.h>

#include <filesystem>

#include "FileNotFoundException.hpp"

namespace coaler {
    class FileParser {
      public:
        static std::vector<RDKit::RWMol *> parse(const std::string &file_path);
    };
}  // namespace coaler
