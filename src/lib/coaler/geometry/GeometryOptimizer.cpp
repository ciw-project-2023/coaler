#include "GeometryOptimizer.hpp"

#include <GraphMol/ChemTransforms/MolFragmenter.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <open3d/Open3D.h>
#include <spdlog/spdlog.h>

#include <iostream>
#include <vector>

namespace coaler {

    GeometryOptimizer::GeometryOptimizer(double max_distance, RDKit::ROMOL_SPTR& core)
        : max_distance_{max_distance}, decomposer_{std::vector<RDKit::ROMOL_SPTR>{core}} {}

    void GeometryOptimizer::save_multi_result_in_vector(multialign::MultiAlignerResult& not_opt_alignment) {
        spdlog::info("Create std::vector< RDKit::ROMOL_SPTR>");
        mol_vec_ = std::vector<RDKit::ROMOL_SPTR>(not_opt_alignment.inputLigands.size());
        pos_id_vec_ = std::vector<multialign::PoseID>(not_opt_alignment.inputLigands.size());
        for (auto id_pair : not_opt_alignment.poseIDsByLigandID) {
            multialign::LigandID ligand_id = id_pair.first;
            multialign::PoseID pose_id = id_pair.second;

            multialign::Ligand ligand = not_opt_alignment.inputLigands.at(ligand_id);

            mol_vec_.at(ligand_id) = boost::make_shared<RDKit::ROMol>(ligand.getMolecule());
            pos_id_vec_.at(ligand_id) = pose_id;
        }
    }

    std::vector<std::vector<open3d::t::geometry::PointCloud>> GeometryOptimizer::create_point_clouds_from_molecules() {
        spdlog::info("Create Point Clouds from Mols");

        // Point Cloud Matrix |rows| == |R_groups + Core|
        std::vector<std::vector<open3d::t::geometry::PointCloud>> point_clouds(
            decomposer_.getRGroupsAsColumns().size(),
            std::vector<open3d::t::geometry::PointCloud>(decomposer_.getRGroupsAsRows().size()));

        size_t group_idx = 0;
        for (auto col : decomposer_.getRGroupsAsColumns()) {
            spdlog::info("R Group {}", col.first);

            // Point Cloud Matrix  |columns| = |mols|
            // std::vector<open3d::t::geometry::PointCloud> point_clouds_r_group(mol_vec_.size());

            size_t mol_idx = 0;
            for (auto mol : col.second) {
                if (RDKit::MolToSmiles(*mol) == "") {
                    spdlog::info("No R-group for mol {} at rgroup {}", mol_idx, group_idx);
                    mol_idx++;
                    continue;
                }

                spdlog::info("Smiles {}", RDKit::MolToSmiles(*mol));

                RDKit::ROMol& mol_o = *mol;
                auto r_group_conf = mol_o.getConformer(pos_id_vec_.at(mol_idx));

                open3d::core::Tensor points
                    = open3d::core::Tensor::Empty({mol->getNumAtoms(), 3}, open3d::core::Float32);

                int i = 0;
                for (auto pos : r_group_conf.getPositions()) {
                    open3d::core::Tensor point = open3d::core::Tensor::Init<double>({pos.x, pos.y, pos.z});
                    points[i] = point;
                    i++;
                }

                open3d::t::geometry::PointCloud cur_cloud(points);
                point_clouds.at(group_idx).at(mol_idx) = cur_cloud;
                mol_idx = mol_idx + 1;
            }
            // point_clouds.at(group_idx) = point_clouds_r_group;
            group_idx = group_idx + 1;
        }

        spdlog::info("Molecule Point Clouds Created");
        return point_clouds;
    }

    void GeometryOptimizer::transform_point_clouds_w_icp(
        std::vector<std::vector<open3d::t::geometry::PointCloud>>& point_clouds) {
        size_t group_idx = 0;
        for (auto col : decomposer_.getRGroupsAsColumns()) {
            spdlog::info("R Group {}", col.first);

            bool target_set = false;
            open3d::t::geometry::PointCloud target_cloud;

            size_t mol_idx = 0;
            for (auto mol : col.second) {
                if (RDKit::MolToSmiles(*mol) == "") {
                    // spdlog::info("No R-group for mol {} at rgroup {}", mol_idx, group_idx);
                    mol_idx++;
                    continue;
                }
                if (!target_set) {
                    spdlog::info("Target cloud set to {} for group {}", mol_idx, group_idx);
                    target_cloud = point_clouds.at(group_idx).at(mol_idx);
                    target_set = true;
                    mol_idx++;
                    continue;
                }

                // spdlog::info("R Group" + std::to_string(col_r));
                // spdlog::info(point_clouds.at(mol).at(col_r).GetPointPositions().ToString());

                auto tranformation = open3d::t::pipelines::registration::ICP(point_clouds.at(group_idx).at(mol_idx),
                                                                             target_cloud, max_distance_);

                point_clouds.at(group_idx).at(mol_idx)
                    = point_clouds.at(group_idx).at(mol_idx).Transform(tranformation.transformation_);

                // spdlog::info("Mol" + std::to_string(mol));
                // spdlog::info(point_clouds.at(mol).at(col_r).GetPointPositions().ToString());

                mol_idx++;
            }
            group_idx++;
        }
        spdlog::info("Molecule Point Clouds Transformed");
    }

    void GeometryOptimizer::set_conformer_pos_to_point_cloud(
        std::vector<std::vector<open3d::t::geometry::PointCloud>>& point_clouds) {
        spdlog::info("Start Setting");

        // rows are molecules and columns are the rgroups
        std::vector<std::vector<RDKit::ROMOL_SPTR>> decomposed(
            decomposer_.getRGroupsAsRows().size(),
            std::vector<RDKit::ROMOL_SPTR>(decomposer_.getRGroupsAsColumns().size()));

        size_t group_idx = 0;
        for (auto col : decomposer_.getRGroupsAsColumns()) {
            //spdlog::info("Group {}", group_idx + 1);
            size_t mol_idx = 0;
            for (auto mol : col.second) {
                if (RDKit::MolToSmiles(*mol) == "") {
                    decomposed.at(mol_idx).at(group_idx) = mol;
                    // spdlog::info("No R-group for mol {} at rgroup {}", mol_idx, group_idx);
                    mol_idx++;
                    continue;
                }

                spdlog::info("Smiles {}", RDKit::MolToSmiles(*mol));

                RDKit::ROMol r_mol = *mol;

                int conformer_count = r_mol.getNumConformers();
                if (conformer_count == 0) {
                    spdlog::info("No Conformere <=> No R Group");
                    decomposed.at(mol_idx).at(group_idx) = boost::make_shared<RDKit::ROMol>(r_mol);
                    mol_idx++;
                    continue;
                }

                // Delete old conformers except one
                auto helper_conf = r_mol.getConformer(pos_id_vec_.at(mol_idx));
                for (int j = 0; j < conformer_count - 1; j++) {
                    auto one_conf = r_mol.getConformer();
                    r_mol.removeConformer(one_conf.getId());
                }
                auto& r_group_conf = r_mol.getConformer();

                // Set to target conformere
                int h_i = 0;
                for (auto& pos : r_group_conf.getPositions()) {
                    pos.x = helper_conf.getPositions()[h_i].x;
                    pos.y = helper_conf.getPositions()[h_i].y;
                    pos.z = helper_conf.getPositions()[h_i].z;
                    h_i++;
                }

                if (!RDKit::MolToSmiles(r_mol).empty() && group_idx != 0) {
                    // spdlog::info("Point Clouds {} {}", point_clouds.size(), point_clouds.at(group_idx).size());
                    spdlog::info("group {} mol {}", group_idx, mol_idx);
                    auto& cur_cloud = point_clouds.at(group_idx).at(mol_idx);
                    //spdlog::info("current cloud get {}", cur_cloud.ToString());
                    auto cur_tensor = cur_cloud.GetPointPositions();
                    // spdlog::info("Current Tensor {} {}", cur_tensor.GetShape(0), cur_tensor.GetShape(1));
                    // spdlog::info("Mol Atoms {}", r_group_conf.getPositions().size());

                    int i = 0;
                    for (auto& pos : r_group_conf.getPositions()) {
                        // spdlog::info("Idx {}", i);
                        pos.x = std::stod(cur_tensor[i][0].ToString());
                        pos.y = std::stod(cur_tensor[i][1].ToString());
                        pos.z = std::stod(cur_tensor[i][2].ToString());
                        i++;
                    }
                }

                decomposed.at(mol_idx).at(group_idx) = boost::make_shared<RDKit::ROMol>(r_mol);

                mol_idx = mol_idx + 1;
            }
            group_idx = group_idx + 1;
        }

//        size_t mol_idx = 0;
//        for (auto row : decomposer_.getRGroupsAsRows()) {
//
//            spdlog::info("MOL");
//
//            std::vector<RDKit::ROMOL_SPTR> decomposed;
//
//            size_t group_idx = 0;
//            for (auto groups : row) {
//                spdlog::info("{} Smiles {}", groups.first, RDKit::MolToSmiles(*groups.second));
//
//                RDKit::ROMol r_mol = *groups.second;
//
//                int conformer_count = r_mol.getNumConformers();
//                if (conformer_count == 0) {
//                    spdlog::info("No Conformere <=> No R Group");
//                    decomposed.emplace_back(boost::make_shared<RDKit::ROMol>(r_mol));
//                    group_idx++;
//                    continue;
//                }
//
//                // Delete old conformers except one
//                auto helper_conf = r_mol.getConformer(pos_id_vec_.at(mol_idx));
//                for (int j = 0; j < conformer_count - 1; j++) {
//                    auto one_conf = r_mol.getConformer();
//                    r_mol.removeConformer(one_conf.getId());
//                }
//                auto& r_group_conf = r_mol.getConformer();
//
//                // Set to target conformere
//                int h_i = 0;
//                for (auto& pos : r_group_conf.getPositions()) {
//                    pos.x = helper_conf.getPositions()[h_i].x;
//                    pos.y = helper_conf.getPositions()[h_i].y;
//                    pos.z = helper_conf.getPositions()[h_i].z;
//                    h_i++;
//                }
//
//                if (!RDKit::MolToSmiles(r_mol).empty()) {
//                    // spdlog::info("Point Clouds {} {}", point_clouds.size(), point_clouds.at(group_idx).size());
//                    // spdlog::info("group {} mol {}", group_idx, mol_idx);
//                    auto& cur_cloud = point_clouds.at(group_idx).at(mol_idx);
//                    // spdlog::info("current cloud get {}", cur_cloud.ToString());
//                    auto cur_tensor = cur_cloud.GetPointPositions();
//                    // spdlog::info("Current Tensor {} {}", cur_tensor.GetShape(0), cur_tensor.GetShape(1));
//                    // spdlog::info("Mol Atoms {}", r_group_conf.getPositions().size());
//
//                    int i = 0;
//                    for (auto& pos : r_group_conf.getPositions()) {
//                        // spdlog::info("Idx {}", i);
//                        pos.x = std::stod(cur_tensor[i][0].ToString());
//                        pos.y = std::stod(cur_tensor[i][1].ToString());
//                        pos.z = std::stod(cur_tensor[i][2].ToString());
//                        i++;
//                    }
//                }
//
//                decomposed.emplace_back(boost::make_shared<RDKit::ROMol>(r_mol));
//                group_idx = group_idx + 1;
//            }
//
//            auto tmp_mol_ = RDKit::molzip(decomposed);
//            RDKit::RWMol mol = *tmp_mol_;
//            spdlog::info("Molzip {}", RDKit::MolToSmiles(mol));
//            geo_opt_ligands_.emplace_back(mol);
//
//            mol_idx = mol_idx + 1;
//        }

        //
        //        size_t group_idx = 0;
        //        for (auto col : decomposer_.getRGroupsAsColumns()) {
        //            spdlog::info("R Group {}", col.first);
        //
        //            size_t mol_idx = 0;
        //            for (auto mol : col.second) {
        //                RDKit::ROMol r_mol = *mol;
        //
        //                int conformer_count = r_mol.getNumConformers();
        //                if(conformer_count == 0){
        //                    spdlog::info("No Conformere <=> No R Group");
        //                    decomposed.at(mol_idx).emplace_back(boost::make_shared<RDKit::ROMol>(r_mol));
        //                    mol_idx++;
        //                    continue;
        //                }
        //
        //                // Delete old conformers except one
        //                auto helper_conf = r_mol.getConformer(pos_id_vec_.at(mol_idx));
        //                for (int j = 0; j < conformer_count - 1; j++) {
        //                    auto one_conf = r_mol.getConformer();
        //                    r_mol.removeConformer(one_conf.getId());
        //                }
        //                auto& r_group_conf = r_mol.getConformer();
        //
        //                // Set to target conformere
        //                int h_i = 0;
        //                for (auto& pos : r_group_conf.getPositions()) {
        //                    pos.x = helper_conf.getPositions()[h_i].x;
        //                    pos.y = helper_conf.getPositions()[h_i].y;
        //                    pos.z = helper_conf.getPositions()[h_i].z;
        //                    h_i++;
        //                }
        //
        //                auto& cur_cloud = point_clouds.at(group_idx).at(mol_idx);
        //                auto& cur_tensor = cur_cloud.GetPointPositions();
        //
        //                if (!RDKit::MolToSmiles(r_mol).empty()) {
        //                    int i = 0;
        //                    for (auto& pos : r_group_conf.getPositions()) {
        //                        // spdlog::info("Before {}", pos.x);
        //                        pos.x = std::stod(cur_tensor[i][0].ToString());
        //                        // spdlog::info("After {}", pos.x);
        //                        pos.y = std::stod(cur_tensor[i][1].ToString());
        //                        pos.z = std::stod(cur_tensor[i][2].ToString());
        //                        i++;
        //                    }
        //                }
        //
        //                decomposed.at(mol_idx).emplace_back(boost::make_shared<RDKit::ROMol>(r_mol));
        //                mol_idx = mol_idx + 1;
        //            }
        //            group_idx = group_idx + 1;
        //        }


        for(int i = 0; i < decomposed.size(); i++){
            auto tmp_mol_ = RDKit::molzip(decomposed.at(i));
            RDKit::RWMol mol = *tmp_mol_;
            spdlog::info("Molzip {}", RDKit::MolToSmiles(mol));

            geo_opt_ligands_.emplace_back(mol);
        }

        spdlog::info("Positions Set");
    }

    void GeometryOptimizer::optimize_alignment_w_icp(multialign::MultiAlignerResult& not_opt_alignment) {
        spdlog::info("Start Geometric Optimization");

        save_multi_result_in_vector(not_opt_alignment);

        for (auto elem : mol_vec_) {
            spdlog::info("Add {}", decomposer_.add(*elem));
        }
        spdlog::info("Process {}", decomposer_.process());

        auto point_clouds = create_point_clouds_from_molecules();

        transform_point_clouds_w_icp(point_clouds);

        set_conformer_pos_to_point_cloud(point_clouds);
    }

    std::vector<RDKit::RWMol> GeometryOptimizer::get_optimized_ligands() { return geo_opt_ligands_; }

    multialign::MultiAlignerResult GeometryOptimizer::get_optimized_alignment() { return geo_opt_alignment_.value(); }
}  // namespace coaler
