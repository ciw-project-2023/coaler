#pragma once

#include <GraphMol/RGroupDecomposition/RGroupDecomp.h>
#include <GraphMol/ROMol.h>
#include <open3d/Open3D.h>

#include "../multialign/MultiAlignerResult.hpp"

namespace coaler {

    class GeometryOptimizer {
      public:
        GeometryOptimizer(double max_distance, RDKit::ROMOL_SPTR& core);

        /**
         *
         * @param not_opt_alignment
         */
        void optimize_alignment_w_icp(multialign::MultiAlignerResult& not_opt_alignment);

        /**
         *
         * @return
         */
        multialign::MultiAlignerResult get_optimized_alignment();

        std::vector<RDKit::RWMol> get_optimized_ligands();

      private:
        void transform_point_clouds_w_icp(std::vector<std::vector<open3d::t::geometry::PointCloud>>& point_clouds);

        void set_conformer_pos_to_point_cloud(std::vector<std::vector<open3d::t::geometry::PointCloud>>&);

        std::vector<std::vector<open3d::t::geometry::PointCloud>> create_point_clouds_from_molecules();

        void save_multi_result_in_vector(multialign::MultiAlignerResult& not_opt_alignment);

        std::optional<multialign::MultiAlignerResult> geo_opt_alignment_{};
        std::vector<RDKit::RWMol> geo_opt_ligands_{};

        std::vector<RDKit::ROMOL_SPTR> mol_vec_{};
        std::vector<multialign::PoseID> pos_id_vec_{};
        // RDKit::RWMol tmp_mol_{};
        // RDKit::Conformer tmp_conf_{};
        RDKit::RGroupDecomposition decomposer_;

        RDKit::RGroupRows rows_{};
        double max_distance_{0};
    };

}  // namespace coaler