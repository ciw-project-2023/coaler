#pragma once

#include <GraphMol/ROMol.h>

#include <unordered_map>

#include "coaler/multialign/Forward.hpp"

namespace coaler::io {
    /**
     * @brief This class is responsible for handling the output.
     */
    class OutputWriter {
      public:
        /**
         * Writes the aligned molecules in an output file.
         * @param file_path
         */
        static void writeSDF(const std::string& filePath, const coaler::multialign::MultiAlignerResult& result);

        /**
         * Writes all conformers of a molecule to a sdf file.
         * @param folder_path A path to the folder store the files in
         * @param mol A molecule whose conformers are to be written to an .sdf file
         */
        static void writeConformersToSDF(const std::string& folderPath, const std::vector<RDKit::ROMOL_SPTR>& mols);
    };

}  // namespace coaler::io
