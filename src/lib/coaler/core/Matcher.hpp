#pragma once

#include <cassert>
#include <queue>
#include <string>
#include <utility>
#include <vector>

#include "GraphMol/FMCS/FMCS.h"
#include "GraphMol/RWMol.h"
#include "GraphMol/Substruct/SubstructMatch.h"
#include "coaler/core/Forward.hpp"
#include "coaler/multialign/models/Forward.hpp"

namespace coaler::core {

    using PairwiseMCSMap = std::unordered_map<multialign::LigandPair,
                                              std::tuple<RDKit::MatchVectType, RDKit::MatchVectType, std::string>,
                                              multialign::LigandPairHash>;

    // NOLINTBEGIN(cppcoreguidelines-pro-type-member-init)
    struct CoreResult {
        RDKit::ROMOL_SPTR core{nullptr};
        RDKit::ROMOL_SPTR ref{nullptr};
        std::unordered_map<int, int> core_to_ref{};
    };
    // NOLINTEND(cppcoreguidelines-pro-type-member-init)

    class Matcher {
      public:
        explicit Matcher(int threads);
        /**
         * calculates the MCS of the molecules
         * @return MCS as ROMol
         */
        std::optional<CoreResult> calculateCoreMcs(RDKit::MOL_SPTR_VECT& mols);

        /**
         * calculates the Murcko scaffold of the molecules
         * @return Murcko Scaffold as ROMol
         */
        std::optional<CoreResult> calculateCoreMurcko(RDKit::MOL_SPTR_VECT& mols);

        /**
         * @return mcs params for very flexibly mcs search
         * i.e. no atom types, bond types, chirality
         */
        static RDKit::MCSParameters getRelaxedMCSParams();

        /**
         * @return mcs params fora strict mcs search
         * i.e. Chirality, Bond order etc.
         */
        static RDKit::MCSParameters getStrictMCSParams();

        static PairwiseMCSMap calcPairwiseMCS(const multialign::LigandVector& mols, bool strict);

      private:
        /**
         * recursive implementation of a murcko pruning of the mcs structure. The function does not delete atoms or
         * bonds from mol but saves them in delAtoms and delBonds to be deleted after function call.
         * @param mol molecule to be pruned
         * @param atomID current atomID of atom looked at
         * @param parentID parent atomID of parent of atom
         * @param visit vector to save visited atoms
         * @param delAtoms vector to save atoms to be deleted
         * @param delBonds vector to save bonds to be deleted
         */
        void murckoPruningRecursive(RDKit::RWMOL_SPTR& mol, int atomID, int parentID, std::vector<bool>& visit,
                                    std::vector<int>& delAtoms, std::vector<std::pair<int, int>>& delBonds,
                                    std::vector<int>& ringAtoms);

        /**
         * recursive implementation to find atoms in delAtoms which are sidechains in the murcko structure and therefore
         * need to stay in the molecule.
         * @param mol molecule to be pruned
         * @param atomID current atomID of atom looked at
         * @param parentID parent atomID of parent of atom
         * @param visit vector to save visited atoms
         * @param ringAtoms atoms which are inside a ring of the molecule  @param mol
         * @param foundRingAtoms atoms found during DFS which are ringatoms of @param mol
         */
        void murckoCheckDelAtoms(RDKit::RWMOL_SPTR& mol, int atomID, int parentID, std::vector<bool>& visit,
                                 std::vector<int>& ringAtoms, std::vector<int>& foundRingAtoms);

        int m_threads;

        [[nodiscard]] RDKit::SubstructMatchParameters getMatchParams() const;

        [[nodiscard]] RDKit::ROMOL_SPTR buildMolConformerForQuery(RDKit::RWMol first, const RDKit::ROMol& query) const;
    };
}  // namespace coaler::core