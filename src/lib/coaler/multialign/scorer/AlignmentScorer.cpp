#include "AlignmentScorer.hpp"

#include <GraphMol/FMCS/FMCS.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/ShapeHelpers/ShapeUtils.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <spdlog/spdlog.h>

namespace coaler::multialign {
    double AlignmentScorer::calcTanimotoShapeSimilarity(const RDKit::ROMol &molA, const RDKit::ROMol &molB,
                                                        unsigned int posIdA, unsigned int posIdB) {
        return 1 - RDKit::MolShapes::tanimotoDistance(molA, molB, posIdA, posIdB);
    }
}  // namespace coaler::multialign