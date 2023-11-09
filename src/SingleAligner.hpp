#pragma once

#include <cstdint> // todo: is necessary when using RDKit libraries

#include <GraphMol/ROMol.h>

namespace coaler {

    /**
     * @brief This class is responsible for alignments of two molecules.
     */
    class SingleAligner {
    public:
        explicit SingleAligner(int core_min_size = 0, int core_max_size = 100);

        /**
         * Computes the alignment of two molecules with a common core structure using the Kabsch' algorithm.
         * If the core structure is not set, then the MCS of the two molecules is calculated.
         * @param mol_a
         * @param mol_b
         * @param core: optional set core structure.
         * @return RMDS score and the core structure of the molecules.
         */
        std::tuple<double, RDKit::ROMOL_SPTR>
        align_molecules_kabsch(RDKit::ROMol mol_a, RDKit::ROMol mol_b, std::optional<RDKit::ROMol> core);


        // TODO: multiple void align_molecules_* functions
    private:
        /**
         * Validates the core structure.
         * @param core
         */
        void validate_core_structure(RDKit::ROMOL_SPTR core) const;

        int core_min_size_{0};
        int core_max_size_{0};
    };

} // coaler
