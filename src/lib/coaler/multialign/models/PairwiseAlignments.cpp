#include "PairwiseAlignments.hpp"

#include <GraphMol/ShapeHelpers/ShapeUtils.h>

#include "Ligand.hpp"

namespace {

    double calc_score(const coaler::multialign::PosePair& key, const std::vector<coaler::multialign::Ligand>& ligands) {
        auto pose1 = key.getFirst();
        auto pose2 = key.getSecond();
        const double distance = RDKit::MolShapes::tanimotoDistance(
            ligands.at(pose1.getLigandId()).getMolecule(), ligands.at(pose2.getLigandId()).getMolecule(),
            static_cast<int>(pose1.getLigandInternalPoseId()), static_cast<int>(pose2.getLigandInternalPoseId()));
        const double similarity = 1 - distance;
        return similarity;
    }
}  // namespace

/*----------------------------------------------------------------------------------------------------------------*/

namespace coaler::multialign {

    double PairwiseAlignments::at(const coaler::multialign::PosePair& key, const LigandVector& ligands, bool store) {
        if (this->count(key) == 1) {
            return this->std::unordered_map<PosePair, double, PosePairHash>::at(key);
        }
        if (!ligands.empty()) {
            if (key == PosePair({0, 29}, {6, 32})) {
            }
            const double score = calc_score(key, ligands);
            if (store) {
                this->emplace(key, score);
            }
            return score;
        }
        throw std::runtime_error("Score not in map and no ligands for calculation provided.");
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    PairwiseAlignments::PairwiseAlignments(PairwiseAlignments& p) : unordered_map(p) {
        for (const auto& elem : p) {
            this->insert(elem);
        }
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    PairwiseAlignments::PairwiseAlignments(const PairwiseAlignments& p) : unordered_map(p) {
        for (const auto& elem : p) {
            this->insert(elem);
        }
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    PairwiseAlignments& PairwiseAlignments::operator=(const PairwiseAlignments& p) {
        this->clear();
        for (const auto& elem : p) {
            this->insert(elem);
        }
        return *this;
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    PairwiseAlignments& PairwiseAlignments::operator=(const std::unordered_map<PosePair, double, PosePairHash>& p) {
        this->clear();
        for (const auto& elem : p) {
            this->insert(elem);
        }
        return *this;
    }

    /*----------------------------------------------------------------------------------------------------------------*/

}  // namespace coaler::multialign