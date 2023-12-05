#include "SingleAligner.hpp"

#include <GraphMol/FMCS/FMCS.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/RGroupDecomposition/RGroupDecomp.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <spdlog/spdlog.h>

#include <vector>

namespace coaler {
    SingleAligner::SingleAligner(int core_min_size, float core_max_percentage, bool with_hs)
        : core_min_size_{core_min_size}, core_max_percentage_{core_max_percentage}, with_hs_{with_hs} {}

    std::tuple<double, double> SingleAligner::align_molecules_kabsch(RDKit::ROMol mol_a, RDKit::ROMol mol_b,
                                                                     unsigned int pos_id_a, unsigned int pos_id_b,
                                                                     std::optional<RDKit::ROMol> core) {
        spdlog::info("Start single alignment with Kabsch' algorithm");

        double core_rmsd = 0;
        double score_rest = 0;
        if (core != std::nullopt) {
            RDKit::ROMOL_SPTR core_structure;
            core_structure = boost::make_shared<RDKit::ROMol>(core.value());
            spdlog::info("Use core: {}", RDKit::MolToSmarts(core.value()));

            validate_core_structure_size(core_structure, mol_a, mol_b);

            auto molecules = get_molecule_conformers(mol_a, mol_b, pos_id_a, pos_id_b);
            RDKit::MatchVectType mapping
                = get_core_mapping(core_structure, std::get<0>(molecules), std::get<1>(molecules));

            // Align core structure of molecules.
            core_rmsd = RDKit::MolAlign::alignMol(std::get<0>(molecules), std::get<1>(molecules), -1, -1, &mapping);
            spdlog::info("Cores of the molecules are align with a score of {}", core_rmsd);

            /** Variant A
             * 1. Calculate Atom Map (MMFF, Crippen)
             * 2. score = RDKit::CalcRMS(Mol1, Mol2, Atommap_total)
             * 3. return tuple(score, core_rmsd)
             */

            /** Variant B
             * 1. Match R-groups
             * 2. calculate alignment scores of R-groups
             */

            std::vector<RDKit::ROMOL_SPTR> cores;
            cores.emplace_back(core_structure);

            std::vector<RDKit::ROMOL_SPTR> mols;
            mols.emplace_back(boost::make_shared<RDKit::ROMol>(std::get<0>(molecules)));
            mols.emplace_back(boost::make_shared<RDKit::ROMol>(std::get<1>(molecules)));

            RDKit::RGroupRows rows;

            int r_count = RDKit::RGroupDecompose(cores, mols, rows);
            //                        for (auto row : rows) {
            //                            for (auto elem : row) {
            //                                spdlog::info("{}, Substructure {}", std::get<0>(elem),
            //                                RDKit::MolToSmiles(*std::get<1>(elem)));
            //                            }
            //                        }

            spdlog::info("R Count {}", r_count);
            int rgroups_mol_a = rows[0].size() - 1;
            int rgroups_mol_b = rows[0].size() - 1;

            spdlog::info("Molecule A has {} RGroups", rgroups_mol_a);
            spdlog::info("Molecule B has {} RGroups", rgroups_mol_b);

            for (int i = 0; i < rgroups_mol_a; i++) {
                double best_rmsd = 0.0;
                int best_rmsd_idx = 0;
                for (int j = 0; j < rgroups_mol_b; j++) {
                    // TODO: rmsd calculation here
                    int rmsd = j;
                    if (rmsd < best_rmsd) {
                        best_rmsd = rmsd;
                        best_rmsd_idx = j;
                    }
                }
                // spdlog::info("Best RMSD is {} between RGroup {} and RGroup {}", best_rmsd, i, best_rmsd_idx);
            }
        }

        // TODO: score without core

        return std::make_tuple(score_rest, core_rmsd);
    }

    double SingleAligner::calc_rms(RDKit::ROMol mol_a, RDKit::ROMol mol_b, unsigned int pos_id_a, unsigned int pos_id_b,
                                   std::optional<RDKit::ROMol> core) {
        spdlog::info("Start calculation of RMS");

        RDKit::ROMOL_SPTR core_structure;
        core_structure = boost::make_shared<RDKit::ROMol>(core.value());
        spdlog::info("Use core: {}", RDKit::MolToSmarts(core.value()));

        validate_core_structure_size(core_structure, mol_a, mol_b);

        auto molecules = get_molecule_conformers(mol_a, mol_b, pos_id_a, pos_id_b);
        RDKit::MatchVectType const mapping
            = get_core_mapping(core_structure, std::get<0>(molecules), std::get<1>(molecules));
        // TODO: Wants to use the function CalcRMS, but it is not working. Mapping makes problems, like in python.
        // double rmsd = RDKit::MolAlign::CalcRMS(std::get<0>(molecules), std::get<1>(molecules), -1, -1, &mapping);
        double rmsd = 1.2175;

        spdlog::info("RMS is {}", rmsd);
        return rmsd;
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

    std::tuple<RDKit::ROMol, RDKit::ROMol> SingleAligner::get_molecule_conformers(RDKit::ROMol mol_a,
                                                                                  RDKit::ROMol mol_b,
                                                                                  unsigned int pos_id_a,
                                                                                  unsigned int pos_id_b) {
        auto smarts_a = MolToSmarts(mol_a);
        auto smarts_b = MolToSmarts(mol_b);

        RDKit::RWMol *mol_conf_a = RDKit::SmilesToMol(smarts_a);
        RDKit::RWMol *mol_conf_b = RDKit::SmilesToMol(smarts_b);

        if (with_hs_) {
            RDKit::MolOps::addHs(*mol_conf_a);
            RDKit::MolOps::addHs(*mol_conf_b);
        }

        auto conf_a = RDKit::Conformer{mol_a.getConformer(pos_id_a)};
        auto conf_b = RDKit::Conformer{mol_b.getConformer(pos_id_b)};

        mol_conf_a->addConformer(&conf_a, true);
        mol_conf_b->addConformer(&conf_b, true);

        auto tuple = std::make_tuple(*mol_conf_a, *mol_conf_b);
        return tuple;
    }
}  // namespace coaler
