//
// Created by chris on 12/7/23.
//

#include "SubstructureAnalyzer.hpp"

#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/new_canon.h>

namespace {
    struct hasDegreeTwo {
        bool operator()(RDKit::Atom* atom) { return atom->getDegree() == 2; }
    };

    /*----------------------------------------------------------------------------------------------------------------*/

    bool allElementsEqual(const std::vector<unsigned> vector) {
        return std::all_of(vector.begin(), vector.end(), [vector](unsigned i) { return i == vector.at(0); });
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    std::vector<unsigned> makeSteadyRankVector(const std::vector<unsigned> originalRanks) {
        std::unordered_map<unsigned, unsigned> originalToSteadyRank;
        unsigned rank = 0;
        for (unsigned orignalRank : originalRanks) {
            if (originalToSteadyRank.count(orignalRank) == 0) {
                originalToSteadyRank.emplace(orignalRank, rank);
                rank++;
            }
        }

        std::vector<unsigned> newRanks;
        for (unsigned originalRank : originalRanks) {
            newRanks.push_back(originalToSteadyRank.at(originalRank));
        }
        return newRanks;
    }
}  // namespace

/*----------------------------------------------------------------------------------------------------------------*/

unsigned coaler::embedder::SubstructureAnalyzer::getNumberOfUniqueSubstructureMatches(const RDKit::ROMol& query,
                                                                                      const RDKit::ROMol& molecule) {
    RDKit::SubstructMatchParameters params;
    params.uniquify = true;
    params.useChirality = true;
    return RDKit::SubstructMatch(molecule, query, params).size();
}

/*----------------------------------------------------------------------------------------------------------------*/

std::pair<unsigned, unsigned> coaler::embedder::SubstructureAnalyzer::getNumberOfRingRotations(const RDKit::ROMol& molecule) {
    // catch mols that have chains
    auto noRotaionsDefault = std::make_pair(1, 0);
    if (!std::all_of(molecule.atoms().begin(), molecule.atoms().end(), hasDegreeTwo())) {
        return noRotaionsDefault;
    }

    unsigned nofAtoms = molecule.getNumAtoms();
    // only even sized rings can be rotation symmetric
    if (nofAtoms % 2 != 0) {
        return noRotaionsDefault;
    }

    std::vector<unsigned> resultCanonIDs;
    RDKit::Canon::rankMolAtoms(molecule, resultCanonIDs, false, true, false);
    std::vector<unsigned> canonIDs = makeSteadyRankVector(resultCanonIDs);

    unsigned max = *std::max_element(canonIDs.begin(), canonIDs.end()) + 1;

    // get number of atoms having rank [index]
    std::vector<unsigned> ranksCount(max);
    for (unsigned rank : canonIDs) {
        ranksCount.at(rank)++;
    }
    /*
     * considers symmetries:
     *  all the same, e.g. c1ccccc1
     *  repeating patterns, e.g. C1NCNCN1 or structures like C1NOSCNOSCNOS1
     */
    unsigned nofRotations = 0;
    switch (max) {
        case 1:
            nofRotations = molecule.getNumAtoms();
            break;
        case 2:
            nofRotations = *std::min_element(ranksCount.begin(), ranksCount.end());
            break;
        default:
            if (allElementsEqual(ranksCount)) {
                return std::make_pair(*ranksCount.begin(), ranksCount.size());
            }
            nofRotations = 1;
    }
    unsigned rotationsStepSize = 0;
    if(nofRotations > 1){
        rotationsStepSize = molecule.getNumAtoms() / nofRotations;
    }
    return std::make_pair(nofRotations,rotationsStepSize);
}
