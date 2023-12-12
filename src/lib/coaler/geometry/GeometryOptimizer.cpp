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

            int i = 0;
            open3d::core::Tensor points = open3d::core::Tensor::Empty({cur_mol.getNumAtoms(), 3}, open3d::core::Float32);
            for (auto pos : cur_conformere.getPositions()) {
                open3d::core::Tensor point = open3d::core::Tensor::Init<double>({pos.x, pos.y, pos.z});
                points[i] = point;
                i++;
            }
            open3d::t::geometry::PointCloud cur_cloud(points);
            point_clouds.emplace_back(cur_cloud);
        }

        spdlog::info("Molecule Point Clouds Created");

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
