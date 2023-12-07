/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <spdlog/spdlog.h>

#include <boost/program_options.hpp>
#include <coaler/embedder/ConformerEmbedder.hpp>
#include <coaler/core/Forward.hpp>
#include <coaler/core/Core.hpp>
#include <coaler/io/Forward.hpp>
#include <coaler/multialign/MultiAligner.hpp>
#include <coaler/multialign/MultiAlignerResult.hpp>
#include <coaler/singlealign/SingleAligner.hpp>
#include <sstream>

namespace opts = boost::program_options;

struct ProgrammOptions {
    std::string input_file_path;
    std::string input_file_type;
    unsigned num_conformers;
    bool dont_add_hydrogens;
    std::string core_type;
};

const std::string help
    = "Usage: aligner [options]\n"
      "Options:\n"
      "  -h, --help\t\t\t\tPrint this help message\n"
      "  -t, --type <type>\t\t\tType of input file (sdf, smiles)\tdefault: smiles\n"
      "  -f, --file <path>\t\t\tPath to input file\n"
      "  -o, --out <path>\t\t\tPath to output file\n"
      "  --core <type>\t\t\tType of core structure used(mcs, murcko)\tdefault: murcko\n"
      "  --conformers <amount>\t\t\tNumber of conformers to generate for each input molecule\tdefault: 10\n"
      "  --dont-add-hydrogens \t\t\tDisabled adding hydrogens to molecules\n";

std::optional<ProgrammOptions> parseArgs(int argc, char* argv[]) {
    ProgrammOptions parsed_options;

    opts::options_description desc("Allowed options");
    desc.add_options()("help,h", "print help message")(
        "type,t", opts::value<std::string>(&parsed_options.input_file_type)->default_value("smiles"),
        "type of input file (sdf, smiles)")(
        "file,f", opts::value<std::string>(&parsed_options.input_file_path)->required(), "path to input file")(
        "core", opts::value<std::string>(&parsed_options.core_type)->default_value("murcko"),
        "type of core structure (mcs, murcko)")(
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

int main(int argc, char* argv[]) {
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

    std::vector<RDKit::RWMol> mols;
    if (opts.input_file_type == "smiles") {
        mols = coaler::io::FileParser::parse(opts.input_file_path);
    } else {
        spdlog::error("provided not supported molecule file format {}", opts.input_file_type);

        return 1;
    }

    spdlog::info("read {} molecules from {} file", mols.size(), opts.input_file_type);

    // generate random core with coordinates. TODO: get coordinates from input

    /*for (auto mol : mols) {
        std::cout << MolToSmiles(mol) << std::endl;
    }
    coaler::core::Core core;
    if (opts.core_type == "mcs") {
        core = coaler::core::Core(mols, coaler::core::CoreType::MCS);
    }
    else {
        core = coaler::core::Core(mols, coaler::core::CoreType::Murcko);
    }*/


    RDKit::RWMol *mol_a = RDKit::SmilesToMol("CCCC");
    RDKit::RWMol *mol_b = RDKit::SmilesToMol("CCCC");
    std::cout << RDKit::MolToSmiles(*mol_a) << std::endl;

    std::vector<RDKit::ROMOL_SPTR> molse;
    molse.emplace_back(boost::make_shared<RDKit::ROMol>(*mol_a));
    molse.emplace_back(boost::make_shared<RDKit::ROMol>(*mol_b));

    RDKit::MCSResult const mcs = RDKit::findMCS(molse);
    RDKit::RWMol coreAsMol = *mcs.QueryMol;

    std::string coreAsString = RDKit::MolToSmiles(coreAsMol);
    std::cout << coreAsString << std::endl;
    std::cout << mcs.SmartsString << std::endl;

    //RDKit::ROMol coreAsMol = core.getCore();
    RDKit::DGeomHelpers::EmbedParameters params;
    RDKit::DGeomHelpers::EmbedMolecule(coreAsMol, params);

    coaler::embedder::CoreAtomMapping coreMapping;

    for (int id = 0; id < coreAsMol.getNumAtoms(); id++) {
        coreMapping.emplace(id, coreAsMol.getConformer(0).getAtomPos(id));
    }

    spdlog::info("embedding {} conformers each into molecules", opts.num_conformers);
    for (RDKit::ROMol& mol : mols) {
        coaler::embedder::ConformerEmbedder conformerEmbedder(coreAsMol, coreMapping);
        if (!conformerEmbedder.embedWithFixedCore(mol, opts.num_conformers)) {
            spdlog::error("Unable to generate conformers. Molecule {} does not match core {}. Aborting.",
                          RDKit::MolToSmiles(mol), RDKit::MolToSmiles(coreAsMol));
            return 1;
        }
    }
    const coaler::SingleAligner singleAligner;
    coaler::multialign::MultiAligner aligner(mols, coreAsMol, singleAligner);
    const coaler::multialign::MultiAlignerResult result = aligner.alignMolecules();

    // write some basic output here to evaluate results

    std::ostringstream oss;
    // takeOwnership must be false for this, as we don't want the SDWriter trying
    // to delete the std::ostringstream.
    bool takeOwnership = false;
    boost::shared_ptr<RDKit::SDWriter> sdf_writer(new RDKit::SDWriter(&oss, takeOwnership));
    if (result.poseIDsByLigandID.size() != result.inputLigands.size()) {
        spdlog::info("only generated an incomplete alignment.");
        return 1;
    }
    for (const auto& entry : result.inputLigands) {
        sdf_writer->write(entry.getMolecule(), result.poseIDsByLigandID.at(entry.getID()));
    }
    std::cout << oss.str() << std::endl;

    spdlog::info("done: exiting");
}
