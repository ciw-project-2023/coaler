/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/MolOps.h>
#include <spdlog/spdlog.h>

#include <boost/program_options.hpp>
#include <coaler/io/Forward.hpp>

namespace opts = boost::program_options;

struct ProgrammOptions {
    std::string input_file_path;
    std::string input_file_type;
    unsigned num_conformers;
    bool dont_add_hydrogens;
};

const std::string help
    = "Usage: aligner [options]\n"
      "Options:\n"
      "  -h, --help\t\t\t\tPrint this help message\n"
      "  -t, --type <type>\t\t\tType of input file (sdf, smiles)\tdefault: smiles\n"
      "  -f, --file <path>\t\t\tPath to input file\n"
      "  -o, --out <path>\t\t\tPath to output file\n"
      "  --conformers <amount>\t\t\tNumber of conformers to generate for each input molecule\tdefault: 10\n"
      "  --dont-add-hydrogens \t\t\tDisabled adding hydrogens to molecules\n";

std::optional<ProgrammOptions> parseArgs(int argc, char *argv[]) {
    ProgrammOptions parsed_options;

    opts::options_description desc("Allowed options");
    desc.add_options()("help,h", "print help message")(
        "type,t", opts::value<std::string>(&parsed_options.input_file_type)->default_value("smiles"),
        "type of input file (sdf, smiles)")(
        "file,f", opts::value<std::string>(&parsed_options.input_file_path)->required(), "path to input file")(
        "out,o", opts::value<std::string>(&parsed_options.input_file_path)->required(), "path to output file")(
        "conformers", opts::value<unsigned>(&parsed_options.num_conformers)->default_value(10))(
        "dont-add-hydrogens", opts::value<bool>(&parsed_options.dont_add_hydrogens)->default_value(false));

    opts::variables_map vm;
    opts::store(opts::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help")) {
        std::cout << help << std::endl;

        return std::nullopt;
    }

    opts::notify(vm);

    return parsed_options;
}

int main(int argc, char *argv[]) {
    std::optional<ProgrammOptions> mOpts;
    try {
        mOpts = parseArgs(argc, argv);
    } catch (...) {
        spdlog::error("failed to parse program arguments");
    }

    if (!mOpts.has_value()) {
        return 0;
    }

    auto opts = mOpts.value();

    std::vector<RDKit::RWMol *> mols;
    if (opts.input_file_type == "smiles") {
        mols = coaler::io::FileParser::parse(opts.input_file_path);
    } else {
        spdlog::error("provided not supported molecule file format {}", opts.input_file_type);

        return 1;
    }

    spdlog::info("read {} molecules from {} file", mols.size(), opts.input_file_type);

    spdlog::info("embedding {} conformers each into molecules", opts.num_conformers);
    for (auto *mol : mols) {
        auto params = RDKit::DGeomHelpers::EmbedParameters{};
        params.randNegEig = false;  // TODO figure out if we need to adjust this and other params here
        RDKit::DGeomHelpers::EmbedMultipleConfs(*mol, opts.num_conformers, params);
    }

    spdlog::info("done: exiting");
}
