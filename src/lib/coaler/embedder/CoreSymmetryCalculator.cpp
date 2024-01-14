//
// Created by chris on 12/6/23.
//

#include "CoreSymmetryCalculator.hpp"

#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
namespace coaler::embedder {

    unsigned CoreSymmetryCalculator::getNofSymmetryAxes(const RDKit::ROMol& mol) {
        std::string const patternString = RDKit::MolToSmarts(mol);
        RDKit::ROMol const patternMol = *RDKit::SmartsToMol(patternString);
        RDKit::SubstructMatchParameters params;
        params.uniquify = false;
        auto matches = RDKit::SubstructMatch(mol, mol, params);
        return matches.size();
    }
}  // namespace coaler::embedder
