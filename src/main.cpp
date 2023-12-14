/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <spdlog/spdlog.h>

#include <boost/program_options.hpp>
#include <coaler/embedder/ConformerEmbedder.hpp>
#include <coaler/geometry/GeometryOptimizer.hpp>
#include <coaler/io/Forward.hpp>
#include <coaler/multialign/MultiAligner.hpp>
#include <coaler/multialign/MultiAlignerResult.hpp>
#include <coaler/singlealign/SingleAligner.hpp>
#include <sstream>

#include "coaler/core/Matcher.hpp"

namespace opts = boost::program_options;
using namespace coaler;

struct ProgrammOptions {
    std::string input_file_path;
    std::string input_file_type;
    unsigned num_conformers;
    std::string out_file;
    int num_threads;
    std::string core_type;
};

const std::string help
    = "Usage: aligner [options]\n"
      "Options:\n"
      "  -h, --help\t\t\t\tPrint this help message\n"
      "  -f, --files <path>\t\t\tPath to input files\n"
      "  -o, --out <path>\t\t\tPath to output files\n"
      "  -j, --threads <amount>\t\t\tNumber of threads to use (default: 1)\n"
      "  --conformers <amount>\t\t\tNumber of conformers to generate for each input molecule (default: 10)\n"
      "  --core <algorithm>\t\t\tAlgorithm to detect core structure (default: mcs, allowed: mcs, murcko)\n";

std::optional<ProgrammOptions> parseArgs(int argc, char* argv[]) {
    ProgrammOptions parsed_options;

    opts::options_description desc("Allowed options");
    desc.add_options()("help,h", "print help message")(
        "file,f", opts::value<std::string>(&parsed_options.input_file_path)->required(), "path to input file")(
        "out,o", opts::value<std::string>(&parsed_options.out_file)->default_value("out.sdf"), "path to output file")(
        "threads,j", opts::value<int>(&parsed_options.num_threads)->default_value(1), "number of threads to use")(
        "core", opts::value<std::string>(&parsed_options.core_type)->default_value("mcs"),
        "algo to detect core structure")("conformers,c",
                                         opts::value<unsigned>(&parsed_options.num_conformers)->default_value(10),
                                         "number of conformers to generate");

    opts::variables_map vm;
    opts::store(opts::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help")) {
        std::cout << help << '\n';

        return std::nullopt;
    }

    opts::notify(vm);

    return parsed_options;
}

int main(int argc, char* argv[]) {
    std::optional<ProgrammOptions> mOpts;
    try {
        mOpts = parseArgs(argc, argv);
    } catch (const boost::program_options::error& e) {
        spdlog::error("failed to parse program arguments: {}", e.what());
    }

    if (!mOpts.has_value()) {
        return 0;
    }

    auto opts = mOpts.value();
    auto mols = io::FileParser::parse(opts.input_file_path);

    std::optional<RDKit::ROMOL_SPTR> core;
    if (opts.core_type == "mcs") {
        core = core::Matcher::calculateCoreMcs(mols);
    } else if (opts.core_type == "murcko") {
        core = core::Matcher::calculateCoreMurcko(mols);
    } else {
        spdlog::error("unknown core calculation algorithm '{}' (allowed values are 'mcs' and 'murcko')",
                      opts.core_type);
        return 1;
    }

    if (!core.has_value()) {
        spdlog::error("failed to calculate core structure");
        return 1;
    }
    // generate random core with coordinates. TODO: get coordinates from input

    const std::string coreSmiles = "C(c1ccccc1)NC";
    RDKit::ROMol* core = RDKit::SmilesToMol(coreSmiles);
    RDKit::DGeomHelpers::EmbedParameters params;
    RDKit::DGeomHelpers::EmbedMolecule(*core, params);


    spdlog::info("embedding {} conformers each into molecules", opts.num_conformers);

    embedder::ConformerEmbedder embedder(core.value(), opts.num_threads);
    for (auto& mol : mols) {
        embedder.embedConformersWithFixedCore(mol, opts.num_conformers);
    }

    spdlog::info("embedding {} conformers each into molecules", opts.num_conformers);
    for (RDKit::ROMol* mol : mols) {
        coaler::embedder::ConformerEmbedder conformerEmbedder(*core, coreMapping);
        if (!conformerEmbedder.embedWithFixedCore(*mol, opts.num_conformers)) {
            spdlog::error("Unable to generate conformers. Molecule {} does not match core {}. Aborting.",
                          RDKit::MolToSmiles(*mol), RDKit::MolToSmiles(*core));
            return 1;
        }
    }
    const coaler::SingleAligner singleAligner;
    coaler::multialign::MultiAligner aligner(mols, *core, singleAligner);
    coaler::multialign::MultiAlignerResult result = aligner.alignMolecules();

    // write some basic output here to evaluate results

    coaler::GeometryOptimizer optimizer(0.005);
    optimizer.optimize_alignment_w_icp(result);
    // coaler::multialign::MultiAlignerResult optimized_result = optimizer.get_optimized_alignment();

    auto ligands = optimizer.get_optimized_ligands();
    spdlog::info("OPTIMIZED LIGANDS {}", ligands.size());

    // TODO: remove

    const std::string file_path = "./optimized_geo.sdf";

    std::ofstream output_file(file_path);
    if (!output_file.is_open()) {
        spdlog::error("Cannot open file: {}", file_path);
        return 1;
    }

    boost::shared_ptr<RDKit::SDWriter> const sdf_writer(new RDKit::SDWriter(&output_file, false));
    for (const auto& entry : ligands) {
        auto conf = entry.getConformer();

        spdlog::info(std::to_string(entry.getNumConformers()));

        for (auto pos : conf.getPositions()) {
            spdlog::info("DONE {}", std::to_string(pos.x));
            spdlog::info("DONE {}", std::to_string(pos.y));
            spdlog::info("DONE {}", std::to_string(pos.z));
        }

        sdf_writer->write(entry);
    }

    //    std::ostringstream oss;
    //    // takeOwnership must be false for this, as we don't want the SDWriter trying
    //    // to delete the std::ostringstream.
    //    bool takeOwnership = false;
    //    boost::shared_ptr<RDKit::SDWriter> sdf_writer(new RDKit::SDWriter(&oss, takeOwnership));
    //    if (result.poseIDsByLigandID.size() != result.inputLigands.size()) {
    //        spdlog::info("only generated an incomplete alignment.");
    //        return 1;
    //    }
    //    for (const auto& entry : result.inputLigands) {
    //        sdf_writer->write(entry.getMolecule(), result.poseIDsByLigandID.at(entry.getID()));
    //    }
    //    std::cout << oss.str() << std::endl;


    spdlog::info("done: exiting");
}
