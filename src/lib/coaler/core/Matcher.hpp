//
// Created by malte on 12/4/23.
//

#pragma once

#include "Forward.hpp"

namespace coaler::core {
    using AtomMap = std::map<int, RDGeom::Point3D>;
    using CoreResult = std::pair<RDKit::ROMOL_SPTR, AtomMap>;

    class Matcher {
      public:
        /**
         * calculates the MCS of the molecules
         * @return MCS as ROMol
         */
        static std::optional<CoreResult> calculateCoreMcs(RDKit::MOL_SPTR_VECT mols, int numOfThreads);

        /**
         * calculates the Murcko scaffold of the molecules
         * @return Murcko Scaffold as ROMol
         */
        static std::optional<CoreResult> calculateCoreMurcko(RDKit::MOL_SPTR_VECT mols, int numOfThreads);

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
        static void murckoPruningRecursive(RDKit::RWMOL_SPTR mol, int atomID, int parentID, std::vector<bool>& visit,
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
        static void murckoCheckDelAtoms(RDKit::RWMOL_SPTR mol, int atomID, int parentID, std::vector<bool>& visit,
                                        std::vector<int>& ringAtoms, std::vector<int>& foundRingAtoms);

    };
}  // namespace coaler::core