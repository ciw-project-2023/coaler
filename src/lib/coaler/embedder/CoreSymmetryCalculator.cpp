#include "CoreSymmetryCalculator.hpp"

#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
namespace coaler::embedder {

    // NOLINTNEXTLINE(readability-convert-member-functions-to-static, misc-unused-parameters) : not our style yet
    [[maybe_unused]] unsigned CoreSymmetryCalculator::getNofSymmetryAxes(const RDKit::ROMol& mol) {
        std::string patternString = RDKit::MolToSmarts(mol);
        RDKit::ROMol patternMol = *RDKit::SmartsToMol(patternString);
        RDKit::SubstructMatchParameters params;
        params.uniquify = false;
        std::vector<RDKit::MatchVectType> matches = RDKit::SubstructMatch(mol, mol, params);
        return matches.size();
    }
}  // namespace coaler::embedder
