#include "OutputWriter.hpp"

#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <spdlog/spdlog.h>

namespace {
    std::string get_formatted_string(unsigned num, unsigned digits) {
        const unsigned numDigits = std::to_string(num).length();
        std::string result;
        while (result.size() < digits - numDigits) {
            result += "0";
        }
        result += std::to_string(num);
        return result;
    }
}  // namespace

/*----------------------------------------------------------------------------------------------------------------*/

namespace coaler::io {
    void OutputWriter::writeSDF(const std::string &filePath, const coaler::multialign::MultiAlignerResult &result) {
        if (result.input_ligands.size() != result.pose_ids_by_ligand_id.size()) {
            throw std::runtime_error(fmt::format("received less output molecules than there was in input: {}/{}",
                                                 result.pose_ids_by_ligand_id.size(), result.input_ligands.size()));
        }
        std::ofstream outputFile(filePath);
        if (!outputFile.is_open()) {
            spdlog::error("Cannot open file: {}", filePath);
            return;
        }

        boost::shared_ptr<RDKit::SDWriter> const sdf_writer(new RDKit::SDWriter(&outputFile, false));
        for (const auto &[ligandId, poseId] : result.pose_ids_by_ligand_id) {
            auto entry = result.input_ligands.at(ligandId).getMolecule();
            entry.setProp("_Score", result.alignment_score);
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

            boost::shared_ptr<RDKit::SDWriter> const sdfWriter(new RDKit::SDWriter(&output_file, false));
            for (unsigned poseId = 0; poseId < mol->getNumConformers(); poseId++) {
                sdfWriter->write(*mol, static_cast<int>(poseId));
            }
        }
    }
}  // namespace coaler::io
