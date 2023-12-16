/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include <boost/program_options.hpp>
#include <coaler/embedder/ConformerEmbedder.hpp>
#include <coaler/io/Forward.hpp>
#include <coaler/multialign/MultiAligner.hpp>
#include <coaler/multialign/MultiAlignerResult.hpp>
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

    std::optional<coaler::core::CoreResult> coreResult;
    if (opts.core_type == "mcs") {
        coreResult = core::Matcher::calculateCoreMcs(mols);
    } else if (opts.core_type == "murcko") {
        coreResult = core::Matcher::calculateCoreMurcko(mols);
    } else {
        spdlog::error("unknown coreResult calculation algorithm '{}' (allowed values are 'mcs' and 'murcko')",
                      opts.core_type);
        return 1;
    }

    if (!coreResult.has_value()) {
        spdlog::error("failed to calculate coreResult structure");
        return 1;
    }

    auto core = coreResult.value();

    spdlog::info("core structure: {}", RDKit::MolToSmarts(*core.first));

    // generate random coreResult with coordinates. TODO: get coordinates from input

    spdlog::info("embedding {} conformers each into molecules", opts.num_conformers);

    embedder::ConformerEmbedder embedder(coreResult->first, coreResult->second, opts.num_threads);
    for (auto& mol : mols) {
        embedder.embedConformersWithFixedCore(mol, opts.num_conformers);
    }

    multialign::MultiAligner aligner(mols);
    auto result = aligner.alignMolecules();

    io::OutputWriter::writeSDF(opts.out_file, result);

    spdlog::info("done: exiting");
}
