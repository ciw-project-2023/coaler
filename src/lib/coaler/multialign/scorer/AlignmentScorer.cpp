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
}  // namespace coaler::multialign