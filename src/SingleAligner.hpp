#pragma once

#include <cstdint> // todo: is necessary when using RDKit libraries

#include <GraphMol/ROMol.h>

namespace ciw {

    /**
     * @brief This class creates an alignment of two molecules.
     */
    class SingleAligner {
    public:
        SingleAligner();

        ~SingleAligner();

        void set_outputfile(std::string);

        void align_molecules_kabsch(RDKit::ROMol mol_a, RDKit::ROMol mol_b, std::optional<RDKit::ROMol> core);

        //void align_molecules_*

    private:
        std::string outputfile;
    };

} // ciw

