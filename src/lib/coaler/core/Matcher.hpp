#pragma once

#include "Forward.hpp"

namespace coaler::core {
    struct CoreResult {
        RDKit::ROMOL_SPTR core{};
        RDKit::ROMOL_SPTR ref{};
        std::unordered_map<int, int> core_to_ref{};
    };

    class Matcher {
      public:
        explicit Matcher(int threads);
        /**
         * calculates the MCS of the molecules
         * @return MCS as ROMol
         */
        std::optional<CoreResult> calculateCoreMcs(const RDKit::MOL_SPTR_VECT& mols);

        /**
         * calculates the Murcko scaffold of the molecules
         * @return Murcko Scaffold as ROMol
         */
        std::optional<CoreResult> calculateCoreMurcko(const RDKit::MOL_SPTR_VECT& mols);

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

        int m_threads;

        [[nodiscard]] RDKit::SubstructMatchParameters getMatchParams() const;

        [[nodiscard]] RDKit::ROMOL_SPTR buildMolConformerForQuery(RDKit::RWMol first, RDKit::ROMol query);
    };
}  // namespace coaler::core