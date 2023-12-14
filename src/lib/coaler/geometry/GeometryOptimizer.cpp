#include "GeometryOptimizer.hpp"

#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <open3d/Open3D.h>
#include <spdlog/spdlog.h>

#include <iostream>
#include <vector>

namespace coaler {

    namespace {
        RDKit::Conformer get_molecule_conformer(RDKit::ROMol mol, unsigned int pos_id) {
            return mol.getConformer(pos_id);
        }
    }  // namespace

    GeometryOptimizer::GeometryOptimizer(double max_distance) : max_distance_{max_distance} {}

    void GeometryOptimizer::optimize_alignment_w_icp(multialign::MultiAlignerResult& not_opt_alignment) {
        spdlog::info("Start Optimizing");

        // Create Point Clouds from Conformer
        std::vector<open3d::t::geometry::PointCloud> point_clouds;
        std::vector<std::tuple<multialign::LigandID, multialign::PoseID>> ligands;
        for (auto id_pair : not_opt_alignment.poseIDsByLigandID) {
            ligands.emplace_back(id_pair);
            auto cur_ligand = not_opt_alignment.inputLigands.at(id_pair.first);
            auto cur_mol = cur_ligand.getMolecule();
            RDKit::Conformer cur_conformere = cur_mol.getConformer(id_pair.second);

            int i = 0;
            open3d::core::Tensor points
                = open3d::core::Tensor::Empty({cur_mol.getNumAtoms(), 3}, open3d::core::Float32);
            for (auto pos : cur_conformere.getPositions()) {
                open3d::core::Tensor point = open3d::core::Tensor::Init<double>({pos.x, pos.y, pos.z});
                points[i] = point;
                i++;
            }
            open3d::t::geometry::PointCloud cur_cloud(points);
            point_clouds.emplace_back(cur_cloud);
        }
        spdlog::info("Molecule Point Clouds Created");

        // Transform Point Clouds using ICP with Target Point Cloud
        open3d::t::geometry::PointCloud target_cloud = point_clouds.at(0);
        for (int idx = 1; idx < point_clouds.size(); idx++) {
            spdlog::info(std::to_string(idx));
            spdlog::info(point_clouds.at(idx).GetPointPositions().ToString());

            auto tranformation
                = open3d::t::pipelines::registration::ICP(point_clouds.at(idx), target_cloud, max_distance_);

            point_clouds.at(idx).Transform(tranformation.transformation_);

            spdlog::info(std::to_string(idx));
            spdlog::info(point_clouds.at(idx).GetPointPositions().ToString());
        }
        spdlog::info("Molecule Point Clouds Transformed");

        // Set Conformere to Point Clouds
        spdlog::info("Size of optimized points {}", point_clouds.size());
        for (int idx = 0; idx < point_clouds.size(); idx++) {
            multialign::Ligand& cur_ligand = not_opt_alignment.inputLigands.at(std::get<0>(ligands.at(idx)));

            auto tmp_mol_ = cur_ligand.getMolecule();
            int conformer_count = tmp_mol_.getNumConformers();
            for (int j = 0; j < conformer_count - 1; j++) {
                auto one_conf = tmp_mol_.getConformer();
                tmp_mol_.removeConformer(one_conf.getId());
            }
            auto& tmp_conf_ = tmp_mol_.getConformer();
            spdlog::info("old conf deleted");

            auto& cur_cloud = point_clouds.at(idx);
            auto& cur_tensor = cur_cloud.GetPointPositions();

            //            auto& positions = tmp_conf_.getPositions();
            //            for(int i = 0; i < positions.size(); i++){
            //                spdlog::info(cur_tensor[i].ToString());
            //                spdlog::info("Before {}", std::to_string(positions[i].x));
            //                positions[i].x = std::stod(cur_tensor[i][0].ToString());// static_cast<double
            //                *>(cur_vector[0].GetDataPtr())[0]; spdlog::info("After {}",
            //                std::to_string(positions[i].x)); positions[i].y =
            //                std::stod(cur_tensor[i][1].ToString());//static_cast<double
            //                *>(cur_vector[1].GetDataPtr())[0]; positions[i].z =
            //                std::stod(cur_tensor[i][2].ToString());//static_cast<double
            //                *>(cur_vector[2].GetDataPtr())[0]; i++;
            //            }
            int i = 0;
            for (auto& pos : tmp_conf_.getPositions()) {
                spdlog::info(cur_tensor[i].ToString());
                spdlog::info("Before {}", std::to_string(pos.x));
                pos.x
                    = std::stod(cur_tensor[i][0].ToString());  // static_cast<double *>(cur_vector[0].GetDataPtr())[0];
                spdlog::info("After {}", std::to_string(pos.x));
                pos.y = std::stod(cur_tensor[i][1].ToString());  // static_cast<double
                                                                 // *>(cur_vector[1].GetDataPtr())[0];
                pos.z = std::stod(cur_tensor[i][2].ToString());  // static_cast<double
                                                                 // *>(cur_vector[2].GetDataPtr())[0];
                i++;
            }

            // tmp_mol_.addConformer(&tmp_conf_);

            geo_opt_ligands_.emplace_back(tmp_mol_);
            for (auto pos : tmp_mol_.getConformer().getPositions()) {
                spdlog::info("GEO DONE {}", std::to_string(pos.x));
                spdlog::info("GEO DONE {}", std::to_string(pos.y));
                spdlog::info("GEO DONE {}", std::to_string(pos.z));
            }

            // cur_ligand.setMolecule(tmp_mol_);
        }

        spdlog::info("DONE");
        // Save in MultiAlignerResult
        // geo_opt_alignment_.emplace(not_opt_alignment.alignmentScore, not_opt_alignment.poseIDsByLigandID,
        // new_ligands);
    }

    std::vector<RDKit::RWMol> GeometryOptimizer::get_optimized_ligands() { return geo_opt_ligands_; }

    multialign::MultiAlignerResult GeometryOptimizer::get_optimized_alignment() { return geo_opt_alignment_.value(); }
}  // namespace coaler
