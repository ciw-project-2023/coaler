#include "AlignmentScorer.hpp"

#include <GraphMol/FMCS/FMCS.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/ShapeHelpers/ShapeUtils.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <spdlog/spdlog.h>

namespace coaler::multialign {
    double AlignmentScorer::calc_tanimoto_shape_similarity(const RDKit::ROMol &mol_a, const RDKit::ROMol &mol_b,
                                                           unsigned int pos_id_a, unsigned int pos_id_b) {
        return 1 - RDKit::MolShapes::tanimotoDistance(mol_a, mol_b, pos_id_a, pos_id_b);
    }

    double AlignmentScorer::calc_alignment_transformation(const RDKit::ROMol &mol_a, const RDKit::ROMol &mol_b,
                                                          unsigned int pos_id_a, unsigned int pos_id_b) {
        RDGeom::Transform3D transformation;
        double rmsd = RDKit::MolAlign::getAlignmentTransform(mol_a, mol_b, transformation, pos_id_a, pos_id_b);
        return rmsd;
    }

    double AlignmentScorer::calc_rmsd(RDKit::ROMol mol_a, RDKit::ROMol mol_b, int pos_id_a, int pos_id_b) {

        auto mols = std::vector<RDKit::ROMOL_SPTR>{RDKit::ROMOL_SPTR(&mol_a), RDKit::ROMOL_SPTR(&mol_b)};
        RDKit::MCSResult const mcs = RDKit::findMCS(mols);

        RDKit::MatchVectType match_vect_a;
        RDKit::SubstructMatch(mol_a, *mcs.QueryMol, match_vect_a);
        RDKit::MatchVectType match_vect_b;
        RDKit::SubstructMatch(mol_b, *mcs.QueryMol, match_vect_b);

        std::vector<RDKit::MatchVectType> mapping;
        mapping.push_back(match_vect_a);
        mapping.push_back(match_vect_b);

        spdlog::info("Calculating RMSD");
        double rmsd = RDKit::MolAlign::CalcRMS(mol_a, mol_b, pos_id_a, pos_id_b, mapping);
        spdlog::info("RMSD: {}", rmsd);
        return rmsd;
    }
}  // namespace coaler::multialign