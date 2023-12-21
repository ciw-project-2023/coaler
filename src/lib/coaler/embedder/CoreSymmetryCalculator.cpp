//
// Created by chris on 12/6/23.
//

#include "CoreSymmetryCalculator.hpp"

#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <spdlog/spdlog.h>
namespace coaler::embedder {

    unsigned CoreSymmetryCalculator::getNofSymmetryAxes(const RDKit::ROMol& mol) {
        std::string patternString = RDKit::MolToSmarts(mol);
        RDKit::ROMol patternMol = *RDKit::SmartsToMol(patternString);
        RDKit::SubstructMatchParameters params;
        params.uniquify = false;
        std::vector<RDKit::MatchVectType> matches = RDKit::SubstructMatch(mol, mol, params);
        return matches.size();
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    CoreAtomMapping CoreSymmetryCalculator::getShiftedMapping(const CoreAtomMapping& map, unsigned int shift) {
        if (shift == 0) {
            return map;
        }
        CoreAtomMapping newMap;
        unsigned idMax = map.size() - 1;
        for (const auto& entry : map) {
            unsigned newId = entry.first + shift;
            if (newId > idMax) {
                newId = newId - idMax - 1;
            }
            newMap.emplace(newId, entry.second);
        }
        return newMap;
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    // TODO add max per match
    std::vector<unsigned> CoreSymmetryCalculator::distributeApproxEvenly(unsigned int nofMatches,
                                                                         unsigned int maxConformers,
                                                                         unsigned int maxPerMatch) {
        std::vector<unsigned> confsForMatch(nofMatches);
        // check if max per match is reached
        if (maxConformers / nofMatches >= maxPerMatch) {
            std::fill(confsForMatch.begin(), confsForMatch.end(), maxPerMatch);
            return confsForMatch;
        }

        // else fill evenly
        unsigned decrementPosition = maxConformers % nofMatches;
        unsigned baseNofConfs = maxConformers / nofMatches;

        std::fill(confsForMatch.begin(), confsForMatch.begin() + decrementPosition, baseNofConfs + 1);

        std::fill(confsForMatch.begin() + decrementPosition, confsForMatch.end(), baseNofConfs);

        return confsForMatch;
    }
}  // namespace coaler::embedder
