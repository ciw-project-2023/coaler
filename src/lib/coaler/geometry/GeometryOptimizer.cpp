#include "GeometryOptimizer.hpp"

#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <open3d/Open3D.h>
#include <spdlog/spdlog.h>

#include <vector>

#include <iostream>

namespace coaler {

    namespace {
        RDKit::Conformer get_molecule_conformer(RDKit::ROMol mol, unsigned int pos_id) {
            return mol.getConformer(pos_id);
        }
    }  // namespace

    GeometryOptimizer::GeometryOptimizer(double max_distance) : max_distance_{max_distance} {}

    void GeometryOptimizer::optimize_alignment_w_icp(multialign::MultiAlignerResult& not_opt_alignment) {
        spdlog::info("Start Optimizing");

        std::vector<open3d::t::geometry::PointCloud> point_clouds;
        for (auto id_pair : not_opt_alignment.poseIDsByLigandID) {
            auto cur_ligand = not_opt_alignment.inputLigands.at(id_pair.first);
            auto cur_mol = cur_ligand.getMolecule();
            RDKit::Conformer cur_conformere = cur_mol.getConformer(id_pair.second);

            open3d::core::Tensor points;
            auto positions = cur_conformere.getPositions();
            for (auto pos : cur_conformere.getPositions()) {
                spdlog::info("Start Point creation");
                //cur_3d_vec.emplace_back(pos.x, pos.y, pos.z);
                open3d::core::Tensor point;
                point.Add(pos.x);
                point.Add(pos.y);
                point.Add(pos.z);
                points.Add(point);
                spdlog::info("Point Created");
            }
            open3d::t::geometry::PointCloud cur_cloud(points);
            point_clouds.emplace_back(cur_cloud);
            spdlog::info("Cloud emplaced");
        }

        spdlog::info("Clouds created");

        open3d::t::geometry::PointCloud target_cloud = point_clouds.at(0);
        for (int idx = 1; idx < point_clouds.size(); idx++) {
            spdlog::info(std::to_string(idx));
            spdlog::info(point_clouds.at(idx).GetPointPositions().ToString());

            auto tranformation
                = open3d::t::pipelines::registration::ICP(point_clouds.at(idx), target_cloud, max_distance_);

            spdlog::info(std::to_string(idx));
            spdlog::info(point_clouds.at(idx).GetPointPositions().ToString());
        }

        // Punktwolke to Mol Conformere

        // Save in MultiAlignerResult

        // set geo_opt_alignment_.emplace(new alignment)
    }

    multialign::MultiAlignerResult GeometryOptimizer::get_optimized_alignment() { return geo_opt_alignment_.value(); }
}  // namespace coaler
