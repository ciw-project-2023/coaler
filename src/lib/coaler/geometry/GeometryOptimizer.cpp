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
            spdlog::info("Start getiing conf");
            auto cur_ligand = not_opt_alignment.inputLigands.at(id_pair.first);
            RDKit::Conformer cur_conformere = cur_ligand.getMolecule().getConformer(id_pair.second); //get_molecule_conformer(cur_ligand.getMolecule(), id_pair.second);
            spdlog::info("conf gettet");

            // point clouds append point cloud of current point moleucle clouds
            //std::vector<Eigen::Vector3d> cur_3d_vec;

            open3d::core::Tensor points;
            spdlog::info("tensor points created");
            const RDGeom::Point3D x = cur_conformere.getAtomPos(0);
            spdlog::info("Conf Positions {}", x.x);
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
