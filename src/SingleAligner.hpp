#pragma once

#include <cstdint> // todo: is necessary when using RDKit libraries

#include <GraphMol/MolAlign/AlignMolecules.h>
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

        void algin_molecules(RDKit::ROMol mol_a, RDKit::ROMol mol_b);

    private:
        std::string outputfile;
    };

} // ciw

