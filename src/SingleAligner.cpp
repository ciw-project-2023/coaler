#include "SingleAligner.hpp"

#include <vector>

#include <spdlog/spdlog.h>

#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/FMCS/FMCS.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>

namespace ciw {
    SingleAligner::SingleAligner(int core_min_size, int core_max_size) : core_min_size_{core_min_size},
                                                                         core_max_size_{core_max_size} {}

    std::tuple<double, RDKit::ROMOL_SPTR>
    SingleAligner::align_molecules_kabsch(RDKit::ROMol mol_a, RDKit::ROMol mol_b, std::optional<RDKit::ROMol> core) {
        /*TODO: Add more conformeres to the molecules with RDKit::DGeomHelpers::EmbedMultipleConfs or Multi-Align
         * has to do these steps in advance, has to be discussed with the group
         * */
        spdlog::info("Start single alignment with Kabsch' algorithm");

        // Add hydrogens to molecules
        RDKit::MolOps::addHs(mol_a);
        RDKit::MolOps::addHs(mol_b);

        // Embed molecules
        RDKit::DGeomHelpers::EmbedMolecule(mol_a);
        RDKit::DGeomHelpers::EmbedMolecule(mol_b);

        RDKit::ROMOL_SPTR core_structure;
        // use mcs if no core structure
        if (!core.has_value()) {
            spdlog::info("No core structure, start calculating MCS");

            // zip mol_a and mol_b into a vector
            std::vector<RDKit::ROMOL_SPTR> mols;
            mols.emplace_back(boost::make_shared<RDKit::ROMol>(mol_a));
            mols.emplace_back(boost::make_shared<RDKit::ROMol>(mol_b));

            // find the MCS
            RDKit::MCSResult res = RDKit::findMCS(mols);
            core_structure = res.QueryMol;
            spdlog::info("MCS: " + res.SmartsString);
        } else {
            core_structure = boost::make_shared<RDKit::ROMol>(core.value());
            spdlog::info("Use core: {}", RDKit::MolToSmarts(core.value()));
        }

        if(core_structure->getNumAtoms() < core_min_size_){
            throw std::runtime_error("Size of core is too small!");
        }
        if(core_structure->getNumAtoms() > core_max_size_){
            throw std::runtime_error("Size of core is too large!");
        }

        // get substructutre match for mol_a and mol_b
        RDKit::MatchVectType match_vect_a;
        RDKit::SubstructMatch(mol_a, *core_structure, match_vect_a);

        RDKit::MatchVectType match_vect_b;
        RDKit::SubstructMatch(mol_b, *core_structure, match_vect_b);

        // zip second value of match_vect_a and match_vect_b into MatchVectType
        RDKit::MatchVectType match_vect;
        for (int i = 0; i < match_vect_a.size(); i++) {
            match_vect.push_back(std::make_pair(match_vect_a[i].second, match_vect_b[i].second));
        }

        // align molecules
        double rmsd = RDKit::MolAlign::alignMol(mol_a, mol_b, -1, -1, &match_vect);

        spdlog::info("Molecules are align with a score of {}", rmsd);

        return std::make_tuple(rmsd, core_structure);
    }

} // ciw
