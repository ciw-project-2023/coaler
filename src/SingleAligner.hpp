#pragma once

#include <cstdint> // todo: is necessary when using RDKit libraries

#include <GraphMol/ROMol.h>

namespace coaler
{

    /**
     * @brief This class creates an alignment of two molecules.
     */
    class SingleAligner
    {
    public:
        explicit SingleAligner(int core_min_size = 0, int core_max_size = 100);

        std::tuple<double, RDKit::ROMOL_SPTR>
        align_molecules_kabsch(RDKit::ROMol mol_a, RDKit::ROMol mol_b, std::optional<RDKit::ROMol> core);

        // void align_molecules_*
    private:
        int core_min_size_{0};
        int core_max_size_{0};
    };

} // ciw
