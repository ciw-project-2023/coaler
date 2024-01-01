/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <omp.h>

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
    unsigned num_start_assemblies;
    std::string core_type;
    bool divideConformersByMatches;
    std::string conformerLogPath;
    bool verbose;
};

const std::string help
    = "Usage: aligner [options]\n"
      "Options:\n"
      "  -h, --help\t\t\t\tPrint this help message\n"
      "  -i, --input-file <path>\t\t\tPath to input files\n"
      "  -o, --out <path>\t\t\tPath to output files\n"
      "  -j, --threads <amount>\t\t\tNumber of threads to use (default: 1)\n"
      "  -v, --verbose\t\t\tActivate verbose logging\n"
      "  --conformers <amount>\t\t\tNumber of conformers per core match to generate for each input molecule (default: "
      "10)\n"
      "  --divide <bool>\t\t\tDivide the number of conformers by the number of times the core is matched in the input "
      "molecule. "
      "Helps against combinatorial explosion if core is small or has high symmetry (default: false)\n"
      "  --assemblies <amount>\t\t\tNumber of starting assemblies (default: 10)\n"
      "  --core <algorithm>\t\t\tAlgorithm to detect core structure (default: mcs, allowed: mcs, murcko)\n"
      "  --confs-log <path>\t\t\tOptional path to folder to store the generated conformers\n";

std::optional<ProgrammOptions> parseArgs(int argc, char* argv[]) {
    ProgrammOptions parsed_options;

    opts::options_description desc("Allowed options");
    desc.add_options()("help,h", "print help message")(
        "input-file,i", opts::value<std::string>(&parsed_options.input_file_path)->required(), "path to input file")(
        "out,o", opts::value<std::string>(&parsed_options.out_file)->default_value("out.sdf"), "path to output file")(
        "threads,j", opts::value<int>(&parsed_options.num_threads)->default_value(1), "number of threads to use")(
        "assemblies, a", opts::value<unsigned>(&parsed_options.num_start_assemblies)->default_value(200),
        "number of starting assemblies to use")(
        "verbose,v", opts::value<bool>(&parsed_options.verbose)->default_value(false), "enable debug logging")(
        "core", opts::value<std::string>(&parsed_options.core_type)->default_value("mcs"),
        "algo to detect core structure")("conformers,c",
                                         opts::value<unsigned>(&parsed_options.num_conformers)->default_value(10),
                                         "number of conformers per core match to generate")(
        "divide,d", opts::value<bool>(&parsed_options.divideConformersByMatches)->default_value(false),
        "divides the number of conformers by the number of times the core is matched")(
        "confs-log", opts::value<std::string>(&parsed_options.conformerLogPath)->default_value("none"));

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
    if (opts.verbose) {
        spdlog::set_level(spdlog::level::debug);
    }

    omp_set_num_threads(opts.num_threads);

    std::ofstream output_file(opts.out_file);
    if (!output_file.is_open()) {
        spdlog::error("Cannot open output file: {}", opts.out_file);
        return 1;
    }

    auto mols = io::FileParser::parse(opts.input_file_path);

    std::optional<coaler::core::CoreResult> coreResult;
    core::Matcher matcher(opts.num_threads);
    spdlog::info("starting core calculation");
    if (opts.core_type == "mcs") {
        coreResult = matcher.calculateCoreMcs(mols);
    } else if (opts.core_type == "murcko") {
        coreResult = matcher.calculateCoreMurcko(mols);
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

    spdlog::info("core structure: {}", RDKit::MolToSmarts(*core.core));

    // generate random coreResult with coordinates. TODO: get coordinates from input

    //const core::PairwiseMCSMap pairwiseStrictMcsMap = matcher.calcPairwiseMCS(mols, true);
    //const core::PairwiseMCSMap pairwiseRelaxedMcsMap = matcher.calcPairwiseMCS(mols, false);

    spdlog::info("embedding {} conformers each into molecules", opts.num_conformers);

    embedder::ConformerEmbedder embedder(core, opts.num_threads, opts.divideConformersByMatches);

#pragma omp parallel for shared(mols, embedder, opts) default(none)
    for (unsigned i = 0; i < mols.size(); i++) {
        embedder.embedConformers(mols.at(i), opts.num_conformers);
    }

    if (opts.conformerLogPath != "none") {
        coaler::io::OutputWriter::writeConformersToSDF(opts.conformerLogPath, mols);
    }

    multialign::MultiAligner aligner(mols,
                                     //pairwiseStrictMcsMap, pairwiseRelaxedMcsMap,
                                     opts.num_start_assemblies, opts.num_threads);
    try {
        auto result = aligner.alignMolecules();
        io::OutputWriter::writeSDF(opts.out_file, result);
    } catch (const std::runtime_error& e)
    {
        spdlog::error(e.what());
        return 5;
    }

    spdlog::info("done: exiting");
}
