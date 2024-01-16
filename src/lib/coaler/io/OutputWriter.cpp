#include "OutputWriter.hpp"

#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <spdlog/spdlog.h>

namespace {
    std::string get_formatted_string(unsigned num, unsigned digits) {
        unsigned num_digits = std::to_string(num).length();
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
    void OutputWriter::writeSDF(const std::string &file_path, const coaler::multialign::MultiAlignerResult &result) {
        if (result.inputLigands.size() != result.poseIDsByLigandID.size()) {
            throw std::runtime_error(fmt::format("received less output molecules than there was in input: {}/{}",
                                                 result.poseIDsByLigandID.size(), result.inputLigands.size()));
        }
        std::ofstream output_file(file_path);
        if (!output_file.is_open()) {
            spdlog::error("Cannot open file: {}", file_path);
            return;
        }

        boost::shared_ptr<RDKit::SDWriter> const sdf_writer(new RDKit::SDWriter(&output_file, false));
        for (const auto &[ligand_id, pose_id] : result.poseIDsByLigandID) {
            auto entry = result.inputLigands.at(ligand_id).getMolecule();
            entry.setProp("_Score", result.alignmentScore);
            sdf_writer->write(entry, pose_id);
        }
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    void OutputWriter::writeConformersToSDF(const std::string &folder_path,
                                            const std::vector<RDKit::ROMOL_SPTR> &mols) {
        std::string file_path = folder_path;
        if (file_path.back() != '/') {
            file_path += "/";
        }
        unsigned magnitude = std::to_string(mols.size()).length();
        for (unsigned molId = 0; molId < mols.size(); molId++) {
            auto mol = mols.at(molId);

            std::string currentFilePath = file_path + "mol_" + get_formatted_string(molId, magnitude) + ".sdf";
            std::ofstream output_file(currentFilePath);

            if (!output_file.is_open()) {
                spdlog::error("Cannot open file: {}", folder_path);
                return;
            }

            boost::shared_ptr<RDKit::SDWriter> const sdf_writer(new RDKit::SDWriter(&output_file, false));
            for (unsigned conformerId = 0; conformerId < mol->getNumConformers(); conformerId++) {
                sdf_writer->write(*mol, conformerId);
            }
        }
    }
}  // namespace coaler::io
