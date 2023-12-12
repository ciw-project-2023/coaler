#pragma once

#include <GraphMol/ROMol.h>

#include "../multialign/MultiAlignerResult.hpp"

namespace coaler {

    class GeometryOptimizer {
      public:
        GeometryOptimizer(double max_distance);

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
        std::optional<multialign::MultiAlignerResult> geo_opt_alignment_{};
        std::vector<RDKit::RWMol> geo_opt_ligands_{};
        //RDKit::RWMol tmp_mol_{};
        //RDKit::Conformer tmp_conf_{};
        double max_distance_{0};
    };

}  // namespace coaler