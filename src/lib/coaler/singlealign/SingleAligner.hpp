/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#pragma once

#include <GraphMol/ROMol.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <cstdint>  // todo: is necessary when using RDKit libraries

namespace coaler {

    /**
     * @brief This class is responsible for alignments of two molecules.
     */
    class SingleAligner {
      public:
        explicit SingleAligner(int core_min_size = 0, float core_max_percentage = 80, bool with_hs = false);

        /**
         * Computes the alignment of two molecules with a common core structure using the Kabsch' algorithm.
         * If the core structure is not set, then the MCS of the two molecules is calculated.
         * @param mol_a
         * @param mol_b
         * @param pos_id_
         * @param pos_id_
         * @param core: optional set core structure.
         * @return RMDS score and the core structure of the molecules.
         */
        double align_molecules_kabsch(RDKit::ROMol mol_a, RDKit::ROMol mol_b, unsigned int pos_id_a,
                                      unsigned int pos_id_b, std::optional<RDKit::ROMol> core);

      private:
        /**
         * Validates the core structure.
         * @param core
         */
        void validate_core_structure_size(RDKit::ROMOL_SPTR core, RDKit::ROMol mol_a, RDKit::ROMol mol_b) const;

        /**
         * Returns the mapping from the subgraph of the core inside molecule a to the subgraph of the core inside
         * molecule b.
         * @param core_structure
         * @param mol_a
         * @param mol_b
         * @return
         */
        RDKit::MatchVectType get_core_mapping(RDKit::ROMOL_SPTR core_structure, RDKit::ROMol mol_a, RDKit::ROMol mol_b);

        /**
         * Returns a tuple of molecules with only one conformer, which is specified with the ids.
         * @param mol_a
         * @param mol_b
         * @param pos_id_a
         * @param pos_id_b
         * @return
         */
        std::tuple<RDKit::ROMol, RDKit::ROMol> get_molecule_conformers(RDKit::ROMol mol_a, RDKit::ROMol mol_b,
                                                                       unsigned int pos_id_a, unsigned int pos_id_b);

        int core_min_size_{0};
        float core_max_percentage_{0};

        bool with_hs_{false};
    };

}  // namespace coaler
