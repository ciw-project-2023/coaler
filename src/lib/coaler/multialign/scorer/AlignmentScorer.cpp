#include "AlignmentScorer.hpp"

#include <GraphMol/FMCS/FMCS.h>
#include <GraphMol/ShapeHelpers/ShapeUtils.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

namespace coaler::multialign {
    double AlignmentScorer::calcTanimotoShapeSimilarity(const RDKit::ROMol &molA, const RDKit::ROMol &molB,
                                                        unsigned int posIdA, unsigned int posIdB) {
        return 1 - RDKit::MolShapes::tanimotoDistance(molA, molB, static_cast<int>(posIdA), static_cast<int>(posIdB));
    }
}  // namespace coaler::multialign