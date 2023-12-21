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

    CoreAtomMapping CoreSymmetryCalculator::getShiftedMapping(const CoreAtomMapping& map,
                                                                                unsigned int shift) {
        CoreAtomMapping newMap;
        unsigned idMax = map.size() - 1;
        for(const auto& entry : map){
            unsigned newId = entry.first + shift;
            if(newId > idMax){
                newId = newId - idMax - 1;
            }
            newMap.emplace(newId, entry.second);
            spdlog::info(newId);
        }
        return newMap;
    }
}  // namespace coaler::embedder
