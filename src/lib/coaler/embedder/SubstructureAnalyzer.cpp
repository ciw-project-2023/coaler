#include "SubstructureAnalyzer.hpp"

#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/new_canon.h>

namespace {
    struct HasDegreeTwo {
        bool operator()(RDKit::Atom* atom) { return atom->getDegree() == 2; }
    };

    /*----------------------------------------------------------------------------------------------------------------*/

    bool AllElementsEqual(const std::vector<unsigned> vector) {
        return std::all_of(vector.begin(), vector.end(), [vector](unsigned i) { return i == vector.at(0); });
    }

    /*----------------------------------------------------------------------------------------------------------------*/
    // NOLINTNEXTLINE(misc-unused-parameters) : clang-tidy no build
    std::vector<unsigned> MakeSteadyRankVector(const std::vector<unsigned> originalRanks) {
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
// NOLINTNEXTLINE(readability-convert-member-functions-to-static) : not our style yet
[[maybe_unused]] unsigned coaler::embedder::SubstructureAnalyzer::getNumberOfUniqueSubstructureMatches(
    // NOLINTNEXTLINE(misc-unused-parameters) : clang-tidy no build
    const RDKit::ROMol& query, const RDKit::ROMol& molecule) {
    RDKit::SubstructMatchParameters params;
    params.uniquify = true;
    params.useChirality = true;
    return RDKit::SubstructMatch(molecule, query, params).size();
}

/*----------------------------------------------------------------------------------------------------------------*/
// NOLINTNEXTLINE(readability-convert-member-functions-to-static) : not our style yet
[[maybe_unused]] unsigned coaler::embedder::SubstructureAnalyzer::getNumberOfRingRotations(
    // NOLINTNEXTLINE(misc-unused-parameters) : clang-tidy no build
    const RDKit::ROMol& molecule) {
    // catch mols that have chains
    if (!std::all_of(molecule.atoms().begin(), molecule.atoms().end(), HasDegreeTwo())) {
        return 1;
    }

    unsigned nofAtoms = molecule.getNumAtoms();
    // only even sized rings can be rotation symmetric
    if (nofAtoms % 2 != 0) {
        return 1;
    }

    std::vector<unsigned> resultCanonIDs;
    RDKit::Canon::rankMolAtoms(molecule, resultCanonIDs, false, true, false);
    std::vector<unsigned> canonIDs = MakeSteadyRankVector(resultCanonIDs);

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
    // NOLINTBEGIN(bugprone-branch-clone) : clang-tidy bug
    switch (max) {
        case 1:
            nofRotations = molecule.getNumAtoms();
            break;
        case 2:
            nofRotations = *std::min_element(ranksCount.begin(), ranksCount.end());
            break;
        default:
            if (AllElementsEqual(ranksCount)) {
                return *ranksCount.begin();
            }
            nofRotations = 1;
    }
    // NOLINTEND(bugprone-branch-clone) : clang-tidy bug

    return nofRotations;
}
