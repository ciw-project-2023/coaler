#include "SingleAligner.hpp"

#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/FMCS/FMCS.h>
#include <vector>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>

namespace ciw {

    SingleAligner::SingleAligner() {
    }

    SingleAligner::~SingleAligner() {
    }

    void SingleAligner::set_outputfile(std::string) {
    }

    void
    SingleAligner::align_molecules_kabsch(RDKit::ROMol mol_a, RDKit::ROMol mol_b, std::optional<RDKit::ROMol> core) {
        /*TODO: Add more conformeres to the molecules with RDKit::DGeomHelpers::EmbedMultipleConfs or Multi-Align
         * has to do these steps in advance, has to be discussed with the group
         * */

        // Add hydrogens to molecules
        RDKit::MolOps::addHs(mol_a);
        RDKit::MolOps::addHs(mol_b);

        // Embed molecules
        RDKit::DGeomHelpers::EmbedMolecule(mol_a);
        RDKit::DGeomHelpers::EmbedMolecule(mol_b);

        RDKit::MCSResult res;
        // use mcs if no core structure
        if (!core.has_value()) {
            // zip mol_a and mol_b into a vector
            std::vector<RDKit::ROMOL_SPTR> mols;
            mols.push_back(RDKit::ROMOL_SPTR(new RDKit::ROMol(mol_a)));
            mols.push_back(RDKit::ROMOL_SPTR(new RDKit::ROMol(mol_b)));

            // find the MCS
            res = RDKit::findMCS(mols);

            // print the MCS
            std::cout << "MCS: " << res.SmartsString << std::endl;
        }
        else {
            //TODO: Check if this really works ... (uff)
            res.SmartsString = core->getProp<std::string>("smarts");
        }
        // Turn core into MCSResult

        // get substructutre match for mol_a and mol_b
        RDKit::MatchVectType match_vect_a;
        RDKit::SubstructMatch(mol_a, *res.QueryMol, match_vect_a);

        RDKit::MatchVectType match_vect_b;
        RDKit::SubstructMatch(mol_b, *res.QueryMol, match_vect_b);

        // zip second value of match_vect_a and match_vect_b into MatchVectType
        RDKit::MatchVectType match_vect;
        for (int i = 0; i < match_vect_a.size(); i++) {
            match_vect.push_back(std::make_pair(match_vect_a[i].second, match_vect_b[i].second));
        }

        // align molecules
        double rmsd = RDKit::MolAlign::alignMol(mol_a, mol_b, -1, -1, &match_vect);

        // print rmsd
        std::cout << "RMSD: " << rmsd << std::endl;
    }

} // ciw