#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <omp.h>

#include <boost/program_options.hpp>
#include <sstream>

#include "coaler/core/Forward.hpp"
#include "coaler/embedder/Forward.hpp"
#include "coaler/io/Forward.hpp"

namespace opts = boost::program_options;
using namespace coaler;

struct ProgrammOptions {
    std::string input_file_path{};
    unsigned num_conformers{};
    std::string out_file{};
    int num_threads{};
    unsigned num_start_assemblies{};
    std::string core_type{};
    bool divide_conformers_by_matches{};
    std::string conformer_log_path{};
    bool verbose{};
    double coarse_optimization_threshold{};
    double fine_optimization_threshold{};
    int optimizer_step_limit{};
};

const std::string HELP
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
      "  --optimizer-coarse_optimization_threshold-threshold <float>\t\t\tTreshold for the "
      "coarse_optimization_threshold "
      "optimization step (default: 1.5)\n"
      "  --optimizer-fine-threshold <float>\t\t\tTreshold for the fine optimization step (default: 0.5)\n"
      "  --optimizer-step-limit <amount> \t\t\tMaximum number of steps for the optimizer (default: 100)\n";

// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays)
std::optional<ProgrammOptions> parse_args(int argc, char* argv[]) {
    ProgrammOptions parsedOptions;

    opts::options_description desc("Allowed options");
    desc.add_options()("help,h", "print help message")(
        "input-file,i", opts::value<std::string>(&parsedOptions.input_file_path)->required(), "path to input file")(
        "out,o", opts::value<std::string>(&parsedOptions.out_file)->default_value("out.sdf"), "path to output file")(
        "threads,j", opts::value<int>(&parsedOptions.num_threads)->default_value(1), "number of threads to use")(
        "assemblies, a", opts::value<unsigned>(&parsedOptions.num_start_assemblies)->default_value(10),
        "number of starting assemblies to use")(
        "verbose,v", opts::value<bool>(&parsedOptions.verbose)->default_value(false), "enable debug logging")(
        "core", opts::value<std::string>(&parsedOptions.core_type)->default_value("mcs"),
        "algo to detect core structure")("conformers,c",
                                         opts::value<unsigned>(&parsedOptions.num_conformers)->default_value(10),
                                         "number of conformers per core match to generate")(
        "divide,d", opts::value<bool>(&parsedOptions.divide_conformers_by_matches)->default_value(false),
        "divides the number of conformers by the number of times the core is matched")(
        "confs-log", opts::value<std::string>(&parsedOptions.conformer_log_path)->default_value("none"))(
        "optimizer-coarse_optimization_threshold-threshold",
        opts::value<double>(&parsedOptions.coarse_optimization_threshold)
            ->default_value(multialign::constants::COARSE_OPTIMIZATION_THRESHOLD))(
        "optimizer-fine-threshold", opts::value<double>(&parsedOptions.fine_optimization_threshold)
                                        ->default_value(multialign::constants::FINE_OPTIMIZATION_THRESHOLD))(
        "optimizer-step-limit", opts::value<int>(&parsedOptions.optimizer_step_limit)
                                    ->default_value(multialign::constants::OPTIMIZER_STEP_LIMIT));

    opts::variables_map vm;
    opts::store(opts::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help")) {
        std::cout << HELP << '\n';

        return std::nullopt;
    }

    opts::notify(vm);

    return parsedOptions;
}

// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays)
int main(int argc, char* argv[]) {
    std::optional<ProgrammOptions> mOpts;
    try {
        mOpts = parse_args(argc, argv);
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
    // NOLINTBEGIN(bugprone-branch-clone)
    if (opts.core_type == "mcs") {
        coreResult = matcher.calculateCoreMcs(mols);
    } else if (opts.core_type == "murcko") {
        coreResult = matcher.calculateCoreMurcko(mols);
    } else {
        spdlog::error("unknown coreResult calculation algorithm '{}' (allowed values are 'mcs' and 'murcko')",
                      opts.core_type);
        return 1;
    }
    // NOLINTEND(bugprone-branch-clone)

    if (!coreResult.has_value()) {
        spdlog::error("failed to calculate coreResult structure");
        return 1;
    }

    auto core = coreResult.value();

    spdlog::info("core structure: {}", RDKit::MolToSmarts(*core.core));

    // generate random coreResult with coordinates. TODO: get coordinates from input

    // const core::PairwiseMCSMap pairwiseStrictMcsMap = matcher.calcPairwiseMCS(mols, true);
    // const core::PairwiseMCSMap pairwiseRelaxedMcsMap = matcher.calcPairwiseMCS(mols, false);

    spdlog::info("Embedding {} conformers for all molecules.", opts.num_conformers);

    embedder::ConformerEmbedder embedder(core, opts.num_threads, opts.divide_conformers_by_matches);

#pragma omp parallel for shared(mols, embedder, opts) default(none)
    for (unsigned i = 0; i < mols.size(); i++) {
        embedder.embedConformers(mols.at(i), opts.num_conformers);
    }

    auto ligands = multialign::LigandVector(mols);
    auto strictMcsMap = coaler::core::Matcher::calcPairwiseMCS(ligands, false);
    auto relaxedMcsMap = coaler::core::Matcher::calcPairwiseMCS(ligands, true);

    const multialign::AssemblyOptimizer optimizer(strictMcsMap, relaxedMcsMap, opts.coarse_optimization_threshold,
                                                  opts.fine_optimization_threshold, opts.optimizer_step_limit,
                                                  opts.num_threads);

    spdlog::info("Finished embedding.");

    if (opts.conformer_log_path != "none") {
        coaler::io::OutputWriter::writeConformersToSDF(opts.conformer_log_path, mols);
    }

    multialign::MultiAligner aligner(mols, optimizer, opts.num_start_assemblies, opts.num_threads);

    const multialign::MultiAlignerResult result = aligner.alignMolecules();
    io::OutputWriter::writeSDF(opts.out_file, result);

    spdlog::info("done: exiting");
}
