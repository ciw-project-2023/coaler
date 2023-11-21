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
        explicit SingleAligner(int core_min_size = 0, float core_max_percentage = 80);

        /**
         * Computes the alignment of two molecules with a common core structure using the Kabsch' algorithm.
         * If the core structure is not set, then the MCS of the two molecules is calculated.
         * @param mol_a
         * @param mol_b
         * @param core: optional set core structure.
         * @return RMDS score and the core structure of the molecules.
         */
        std::tuple<double, RDKit::ROMOL_SPTR> align_molecules_kabsch(RDKit::ROMol mol_a, RDKit::ROMol mol_b, std::optional<RDKit::ROMol> core);

        // TODO: multiple void align_molecules_* functions
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

        int core_min_size_{0};
        float core_max_percentage_{0};
    };

}  // namespace coaler
