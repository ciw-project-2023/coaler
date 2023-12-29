#include "AlignmentScorer.hpp"

#include <GraphMol/FMCS/FMCS.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/ShapeHelpers/ShapeUtils.h>
#include <GraphMol/SmilesParse/SmilesParse.h>


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
}  // namespace coaler