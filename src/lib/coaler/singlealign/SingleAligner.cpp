/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#include "SingleAligner.hpp"

#include <GraphMol/FMCS/FMCS.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/ShapeHelpers/ShapeUtils.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <spdlog/spdlog.h>

#include <vector>

namespace coaler {
    SingleAligner::SingleAligner(int core_min_size, float core_max_percentage, bool with_hs)
        : core_min_size_{core_min_size}, core_max_percentage_{core_max_percentage}, with_hs_{with_hs} {
        if (with_hs_) {
            spdlog::info("SingleAligner is initiated, minimum core size set to {} and does consider H-Atoms ",
                         core_min_size_);
        } else {
            spdlog::info("SingleAligner is initiated, minimum core size set to {} and does not consider H-Atoms ",
                         core_min_size_);
        }
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    double SingleAligner::align_molecules_kabsch(const RDKit::ROMol& mol_a, const RDKit::ROMol& mol_b,
                                                 unsigned int pos_id_a, unsigned int pos_id_b, RDKit::ROMol core) {
        RDKit::ROMOL_SPTR core_structure;
        core_structure = boost::make_shared<RDKit::ROMol>(core);

        validate_core_structure_size(core_structure, mol_a, mol_b);

        auto molecules = get_molecule_conformers(mol_a, mol_b, pos_id_a, pos_id_b);
        RDKit::MatchVectType mapping = get_core_mapping(core_structure, std::get<0>(molecules), std::get<1>(molecules));

        return RDKit::MolAlign::alignMol(std::get<0>(molecules), std::get<1>(molecules), -1, -1, &mapping);
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    double SingleAligner::calculate_tanimoto_shape_similarity(const RDKit::ROMol& mol_a, const RDKit::ROMol& mol_b,
                                                              unsigned int pos_id_a, unsigned int pos_id_b) {
        // auto molecules = get_molecule_conformers(mol_a, mol_b, pos_id_a, pos_id_b);
        // return 1 - RDKit::MolShapes::tanimotoDistance(std::get<0>(molecules), std::get<1>(molecules));
        return 1 - RDKit::MolShapes::tanimotoDistance(mol_a, mol_b, pos_id_a, pos_id_b);
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    void SingleAligner::validate_core_structure_size(RDKit::ROMOL_SPTR core, const RDKit::ROMol& mol_a,
                                                     const RDKit::ROMol& mol_b) const {
        if (core->getNumAtoms() < core_min_size_) {
            spdlog::error("Size of core is too small!");
            throw std::runtime_error("Size of core is too small!");
        }

        float core_percentage_to_mol_a = (core->getNumAtoms() * 100) / mol_a.getNumAtoms();
        float core_percentage_to_mol_b = (core->getNumAtoms() * 100) / mol_b.getNumAtoms();
        if (core_percentage_to_mol_a > core_max_percentage_ || core_percentage_to_mol_b > core_max_percentage_) {
            spdlog::warn("Size of core is near the size of the input molecules!");
        }
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    RDKit::MatchVectType SingleAligner::get_core_mapping(RDKit::ROMOL_SPTR core_structure, const RDKit::ROMol& mol_a,
                                                         const RDKit::ROMol& mol_b) {
        // Find core inside molecules.
        RDKit::MatchVectType match_vect_a;
        RDKit::SubstructMatch(mol_a, *core_structure, match_vect_a);

        RDKit::MatchVectType match_vect_b;
        RDKit::SubstructMatch(mol_b, *core_structure, match_vect_b);

        if (match_vect_a.size() != core_structure->getNumAtoms()
            || match_vect_b.size() != core_structure->getNumAtoms()) {
            spdlog::error("Core is not a common core structure of molecule a and molecule b!");
            throw std::runtime_error("Core is not a common core structure of molecule a and molecule b!");
        }

        RDKit::MatchVectType mapping;
        for (int i = 0; i < match_vect_a.size(); i++) {
            mapping.push_back(std::make_pair(match_vect_a[i].second, match_vect_b[i].second));
        }
        return mapping;
    }
}  // namespace coaler
