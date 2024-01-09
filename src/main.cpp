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
    std::string inputFilePath;
    unsigned numConformers;
    std::string outFile;
    int numThreads;
    unsigned numStartAssemblies;
    std::string coreType;
    bool divideConformersByMatches;
    std::string conformerLogPath;
    bool verbose;
    double coarseOptimizationThreshold;
    double fineOptimizationThreshold;
    int optimizerStepLimit;
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
      "  --confs-log <path>\t\t\tOptional path to folder to store the generated conformers\n"
      "  --optimizer-coarseOptimizationThreshold-threshold <float>\t\t\tTreshold for the coarseOptimizationThreshold "
      "optimization step (default: 1.5)\n"
      "  --optimizer-fine-threshold <float>\t\t\tTreshold for the fine optimization step (default: 0.5)\n"
      "  --optimizer-step-limit <amount> \t\t\tMaximum number of steps for the optimizer (default: 100)\n";

std::optional<ProgrammOptions> parseArgs(int argc, char* argv[]) {
    ProgrammOptions parsed_options;

    opts::options_description desc("Allowed options");
    desc.add_options()("help,h", "print help message")(
        "input-file,i", opts::value<std::string>(&parsed_options.inputFilePath)->required(), "path to input file")(
        "out,o", opts::value<std::string>(&parsed_options.outFile)->default_value("out.sdf"), "path to output file")(
        "threads,j", opts::value<int>(&parsed_options.numThreads)->default_value(1), "number of threads to use")(
        "assemblies, a", opts::value<unsigned>(&parsed_options.numStartAssemblies)->default_value(10),
        "number of starting assemblies to use")(
        "verbose,v", opts::value<bool>(&parsed_options.verbose)->default_value(false), "enable debug logging")(
        "core", opts::value<std::string>(&parsed_options.coreType)->default_value("mcs"),
        "algo to detect core structure")("conformers,c",
                                         opts::value<unsigned>(&parsed_options.numConformers)->default_value(10),
                                         "number of conformers per core match to generate")(
        "divide,d", opts::value<bool>(&parsed_options.divideConformersByMatches)->default_value(false),
        "divides the number of conformers by the number of times the core is matched")(
        "confs-log", opts::value<std::string>(&parsed_options.conformerLogPath)->default_value("none"))(
        "optimizer-coarseOptimizationThreshold-threshold",
        opts::value<double>(&parsed_options.coarseOptimizationThreshold)
            ->default_value(multialign::Constants::COARSE_OPTIMIZATION_THRESHOLD))(
        "optimizer-fine-threshold", opts::value<double>(&parsed_options.fineOptimizationThreshold)
                                        ->default_value(multialign::Constants::FINE_OPTIMIZATION_THRESHOLD))(
        "optimizer-step-limit", opts::value<int>(&parsed_options.optimizerStepLimit)
                                    ->default_value(multialign::Constants::OPTIMIZER_STEP_LIMIT));

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

    omp_set_num_threads(opts.numThreads);

    std::ofstream output_file(opts.outFile);
    if (!output_file.is_open()) {
        spdlog::error("Cannot open output file: {}", opts.outFile);
        return 1;
    }

    auto mols = io::FileParser::parse(opts.inputFilePath);

    std::optional<coaler::core::CoreResult> coreResult;
    core::Matcher matcher(opts.numThreads);
    spdlog::info("starting core calculation");
    if (opts.coreType == "mcs") {
        coreResult = matcher.calculateCoreMcs(mols);
    } else if (opts.coreType == "murcko") {
        coreResult = matcher.calculateCoreMurcko(mols);
    } else {
        spdlog::error("unknown coreResult calculation algorithm '{}' (allowed values are 'mcs' and 'murcko')",
                      opts.coreType);
        return 1;
    }

    if (!coreResult.has_value()) {
        spdlog::error("failed to calculate coreResult structure");
        return 1;
    }

    auto core = coreResult.value();

    spdlog::info("core structure: {}", RDKit::MolToSmarts(*core.core));

    // generate random coreResult with coordinates. TODO: get coordinates from input

    // const core::PairwiseMCSMap pairwiseStrictMcsMap = matcher.calcPairwiseMCS(mols, true);
    // const core::PairwiseMCSMap pairwiseRelaxedMcsMap = matcher.calcPairwiseMCS(mols, false);

    spdlog::info("Embedding {} conformers for all molecules.", opts.numConformers);

    embedder::ConformerEmbedder embedder(core, opts.numThreads, opts.divideConformersByMatches);

#pragma omp parallel for shared(mols, embedder, opts) default(none)
    for (unsigned i = 0; i < mols.size(); i++) {
        embedder.embedConformers(mols.at(i), opts.numConformers);
    }

    auto ligands = multialign::LigandVector(mols);
    auto strictMcsMap = coaler::core::Matcher::calcPairwiseMCS(ligands, false);
    auto relaxedMcsMap = coaler::core::Matcher::calcPairwiseMCS(ligands, true);

    const multialign::AssemblyOptimizer optimizer(strictMcsMap, relaxedMcsMap, opts.coarseOptimizationThreshold,
                                                  opts.fineOptimizationThreshold, opts.optimizerStepLimit,
                                                  opts.numThreads);

    spdlog::info("Finished embedding.");

    if (opts.conformerLogPath != "none") {
        coaler::io::OutputWriter::writeConformersToSDF(opts.conformerLogPath, mols);
    }

    multialign::MultiAligner aligner(mols, optimizer, opts.numStartAssemblies, opts.numThreads);

    const multialign::MultiAlignerResult result = aligner.alignMolecules();
    io::OutputWriter::writeSDF(opts.outFile, result);

    spdlog::info("done: exiting");
}
