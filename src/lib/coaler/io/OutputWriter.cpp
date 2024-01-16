#include "OutputWriter.hpp"

#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <spdlog/spdlog.h>

namespace {
    std::string get_formatted_string(unsigned num, unsigned digits) {
        const unsigned num_digits = std::to_string(num).length();
        std::string result;
        while (result.size() < digits - num_digits) {
            result += "0";
        }
        result += std::to_string(num);
        return result;
    }
}  // namespace

/*----------------------------------------------------------------------------------------------------------------*/

namespace coaler::io {
    void OutputWriter::writeSDF(const std::string &filePath, const coaler::multialign::MultiAlignerResult &result) {
        if (result.inputLigands.size() != result.poseIDsByLigandID.size()) {
            throw std::runtime_error(fmt::format("received less output molecules than there was in input: {}/{}",
                                                 result.poseIDsByLigandID.size(), result.inputLigands.size()));
        }
        std::ofstream output_file(filePath);
        if (!output_file.is_open()) {
            spdlog::error("Cannot open file: {}", filePath);
            return;
        }

        boost::shared_ptr<RDKit::SDWriter> const sdf_writer(new RDKit::SDWriter(&output_file, false));
        for (const auto &[ligandId, poseId] : result.poseIDsByLigandID) {
            auto entry = result.inputLigands.at(ligandId).getMolecule();
            entry.setProp("_Score", result.alignmentScore);
            sdf_writer->write(entry, static_cast<int>(poseId));
        }
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    void OutputWriter::writeConformersToSDF(const std::string &folderPath, const std::vector<RDKit::ROMOL_SPTR> &mols) {
        std::string filePath = folderPath;
        if (filePath.back() != '/') {
            filePath += "/";
        }
        const unsigned magnitude = std::to_string(mols.size()).length();
        for (unsigned molId = 0; molId < mols.size(); molId++) {
            const auto &mol = mols.at(molId);

            const std::string currentFilePath = filePath + "mol_" + get_formatted_string(molId, magnitude) + ".sdf";
            std::ofstream output_file(currentFilePath);

            if (!output_file.is_open()) {
                spdlog::error("Cannot open file: {}", folderPath);
                return;
            }

            boost::shared_ptr<RDKit::SDWriter> const sdf_writer(new RDKit::SDWriter(&output_file, false));
            for (unsigned poseId = 0; poseId < mol->getNumConformers(); poseId++) {
                sdf_writer->write(*mol, static_cast<int>(poseId));
            }
        }
    }
}  // namespace coaler::io
