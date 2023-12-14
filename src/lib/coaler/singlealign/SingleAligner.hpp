/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#pragma once

#include <GraphMol/ROMol.h>
#include <GraphMol/Substruct/SubstructMatch.h>

namespace coaler {

    /**
     * @brief This class is responsible for alignments of two molecules.
     */
    class SingleAligner {
      public:
        explicit SingleAligner(int core_min_size = 0, float core_max_percentage = 80, bool with_hs = false);

        /**
         * Computes the RMSD of the core inside the molecule using the Kabsch' algorithm.
         * @param mol_a
         * @param mol_b
         * @param pos_id_a: Conforemere ID for molecules A
         * @param pos_id_b: Conforemere ID for molecules B
         * @param core
         * @return RMDS score of the core inside the molecules.
         */
        double align_molecules_kabsch(const RDKit::ROMol& mol_a, const RDKit::ROMol& mol_b, unsigned int pos_id_a,
                                      unsigned int pos_id_b, RDKit::ROMol core);

        /**
         * Computes the tanimoto shape similarity between two molecules.
         * @param mol_a
         * @param mol_b
         * @param pos_id_a: Conforemere ID for molecules A
         * @param pos_id_b: Conforemere ID for molecules B
         * @return Tanimoto shape similarity of molecules
         */
        double calculate_tanimoto_shape_similarity(const RDKit::ROMol& mol_a, const RDKit::ROMol& mol_b,
                                                   unsigned int pos_id_a, unsigned int pos_id_b);

      private:
        /**
         * Validates the core structure.
         * @param core
         */
        void validate_core_structure_size(RDKit::ROMOL_SPTR core, const RDKit::ROMol& mol_a,
                                          const RDKit::ROMol& mol_b) const;

        /**
         * Returns the mapping from the subgraph of the core inside molecule a to the subgraph of the core inside
         * molecule b.
         * @param core_structure
         * @param mol_a
         * @param mol_b
         * @return
         */
        RDKit::MatchVectType get_core_mapping(RDKit::ROMOL_SPTR core_structure, const RDKit::ROMol& mol_a,
                                              const RDKit::ROMol& mol_b);

        int core_min_size_{0};
        float core_max_percentage_{0};

        bool with_hs_{false};
    };

}  // namespace coaler
