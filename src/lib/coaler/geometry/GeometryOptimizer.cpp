#include "GeometryOptimizer.hpp"

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <open3d/Open3D.h>
#include <spdlog/spdlog.h>

#include <iostream>
#include <vector>

namespace coaler {

    namespace {
        std::vector<open3d::t::geometry::PointCloud> create_point_clouds_from_molecules(
            multialign::MultiAlignerResult& not_opt_alignment) {
            std::vector<open3d::t::geometry::PointCloud> point_clouds;

            for (auto id_pair : not_opt_alignment.poseIDsByLigandID) {
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

            return point_clouds;
        }
    }  // namespace

    GeometryOptimizer::GeometryOptimizer(double max_distance) : max_distance_{max_distance} {}

    void GeometryOptimizer::transform_point_clouds_w_icp(std::vector<open3d::t::geometry::PointCloud>& point_clouds) {
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
    }

    void GeometryOptimizer::set_conformer_pos_to_point_cloud(multialign::Ligand& cur_ligand,
                                                             std::vector<open3d::t::geometry::PointCloud>& point_clouds,
                                                             size_t idx) {
        auto tmp_mol_ = cur_ligand.getMolecule();
        int conformer_count = tmp_mol_.getNumConformers();
        for (int j = 0; j < conformer_count - 1; j++) {
            auto one_conf = tmp_mol_.getConformer();
            tmp_mol_.removeConformer(one_conf.getId());
        }
        auto& tmp_conf_ = tmp_mol_.getConformer();

        auto& cur_cloud = point_clouds.at(idx);
        auto& cur_tensor = cur_cloud.GetPointPositions();

        int i = 0;
        for (auto& pos : tmp_conf_.getPositions()) {
            pos.x = std::stod(cur_tensor[i][0].ToString());
            pos.y = std::stod(cur_tensor[i][1].ToString());
            pos.z = std::stod(cur_tensor[i][2].ToString());
            i++;
        }

        geo_opt_ligands_.emplace_back(tmp_mol_);
    }

    void GeometryOptimizer::optimize_alignment_w_icp(multialign::MultiAlignerResult& not_opt_alignment) {
        spdlog::info("Start Optimizing");

        std::vector<std::tuple<multialign::LigandID, multialign::PoseID>> ligands;
        for (auto id_pair : not_opt_alignment.poseIDsByLigandID) {
            ligands.emplace_back(id_pair);
        }

        auto point_clouds = create_point_clouds_from_molecules(not_opt_alignment);

        transform_point_clouds_w_icp(point_clouds);

        // Set Conformere to Point Clouds
        for (int idx = 0; idx < point_clouds.size(); idx++) {
            multialign::Ligand& cur_ligand = not_opt_alignment.inputLigands.at(std::get<0>(ligands.at(idx)));
            set_conformer_pos_to_point_cloud(cur_ligand, point_clouds, idx);
        }
    }

    std::vector<RDKit::RWMol> GeometryOptimizer::get_optimized_ligands() { return geo_opt_ligands_; }

    multialign::MultiAlignerResult GeometryOptimizer::get_optimized_alignment() { return geo_opt_alignment_.value(); }
}  // namespace coaler
