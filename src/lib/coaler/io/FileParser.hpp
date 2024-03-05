#include <GraphMol/SmilesParse/SmilesParse.h>
#include <spdlog/spdlog.h>

#include <filesystem>

#include "FileNotFoundException.hpp"

namespace coaler::io {
    class FileParser {
        /*!
         * @brief Parses a file and returns a vector of molecules
         *
         * @param filePath Path to the file
         * @return RDKit::MOL_SPTR_VECT Vector of molecules
         */
      public:
        /*!
         * @brief Parses a file and returns a vector of molecules
         *
         * @param filePath Path to the file
         * @return RDKit::MOL_SPTR_VECT Vector of molecules
         */
        static RDKit::MOL_SPTR_VECT parse(const std::string& filePath);

      private:
        /*!
         * @brief Checks if the input molecules are valid
         *
         * @param mols Vector of molecules
         * @param filePath Path to the file
         * @return RDKit::MOL_SPTR_VECT Vector of molecules
         */
        static RDKit::MOL_SPTR_VECT checkInputMolecules(const RDKit::MOL_SPTR_VECT& mols, const std::string& filePath);
    };
}  // namespace coaler::io
