#include "SingleAligner.hpp"

#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/FMCS/FMCS.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <spdlog/spdlog.h>

#include <vector>

namespace coaler {
    SingleAligner::SingleAligner(int core_min_size, float core_max_percentage)
        : core_min_size_{core_min_size}, core_max_percentage_{core_max_percentage} {}

    std::tuple<double, RDKit::ROMOL_SPTR> SingleAligner::align_molecules_kabsch(RDKit::ROMol mol_a, RDKit::ROMol mol_b,
                                                                                std::optional<RDKit::ROMol> core) {
        /*TODO: Add more conformeres to the molecules with RDKit::DGeomHelpers::EmbedMultipleConfs or
         * Multi-Align has to do these steps in advance, has to be discussed with the group
         * */
        spdlog::info("Start single alignment with Kabsch' algorithm");

        RDKit::MolOps::addHs(mol_a);
        RDKit::MolOps::addHs(mol_b);

        RDKit::DGeomHelpers::EmbedMolecule(mol_a);
        RDKit::DGeomHelpers::EmbedMolecule(mol_b);

        RDKit::ROMOL_SPTR core_structure;
        // calculate mcs if no core structure
        if (!core.has_value()) {
            spdlog::info("No core structure, start calculating MCS");

            // zip mol_a and mol_b into a vector
            std::vector<RDKit::ROMOL_SPTR> mols;
            mols.emplace_back(boost::make_shared<RDKit::ROMol>(mol_a));
            mols.emplace_back(boost::make_shared<RDKit::ROMol>(mol_b));

            RDKit::MCSResult res = RDKit::findMCS(mols);
            core_structure = res.QueryMol;
            spdlog::info("MCS: " + res.SmartsString);
        } else {
            core_structure = boost::make_shared<RDKit::ROMol>(core.value());
            spdlog::info("Use core: {}", RDKit::MolToSmarts(core.value()));
        }

        validate_core_structure_size(core_structure, mol_a, mol_b);

        RDKit::MatchVectType mapping = get_core_mapping(core_structure, mol_a, mol_b);
        double rmsd = RDKit::MolAlign::alignMol(mol_a, mol_b, -1, -1, &mapping);

        spdlog::info("Molecules are align with a score of {}", rmsd);
        return std::make_tuple(rmsd, core_structure);
    }

    void SingleAligner::validate_core_structure_size(RDKit::ROMOL_SPTR core, RDKit::ROMol mol_a,
                                                     RDKit::ROMol mol_b) const {
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

    RDKit::MatchVectType SingleAligner::get_core_mapping(RDKit::ROMOL_SPTR core_structure, RDKit::ROMol mol_a,
                                                         RDKit::ROMol mol_b) {
        // find core inside molecules
        RDKit::MatchVectType match_vect_a;
        RDKit::SubstructMatch(mol_a, *core_structure, match_vect_a);

        RDKit::MatchVectType match_vect_b;
        RDKit::SubstructMatch(mol_b, *core_structure, match_vect_b);

        if (match_vect_a.size() != match_vect_b.size()) {
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
