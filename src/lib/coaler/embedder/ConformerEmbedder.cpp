/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#include "ConformerEmbedder.hpp"

#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/FMCS/FMCS.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <spdlog/spdlog.h>

#include <utility>

const unsigned seed = 42;
const float forceTol = 0.0135;

namespace {

    /*----------------------------------------------------------------------------------------------------------------*/

    RDKit::DGeomHelpers::EmbedParameters get_embed_params_for_optimizer_generation() {
        RDKit::DGeomHelpers::EmbedParameters params;
        params = RDKit::DGeomHelpers::srETKDGv3;
        params.optimizerForceTol = forceTol;
        params.randomSeed = seed;
        params.useRandomCoords = true;
        params.numThreads = 1;
        params.clearConfs = false;
        return params;
    }
}  // namespace

namespace coaler::embedder {
    ConformerEmbedder::ConformerEmbedder(core::CoreResult result, const int threads,
                                         const bool divideConformersByMatches)
        : m_core(std::move(result)), m_threads(threads), m_divideConformersByMatches(divideConformersByMatches) {}

    void ConformerEmbedder::embedConformers(const RDKit::ROMOL_SPTR &mol, multialign::LigandID id, unsigned numConfs) {
        auto relaxedMcs = m_core.pairwiseMCSResult.relaxedMcsMap;
        if ( id == m_core.refID) {
           return;
        }

        multialign::LigandPair key(id, m_core.refID);
        auto [molMatch, refMatch, _] = relaxedMcs.at(key);

        spdlog::debug("number of Core Matches: {}", molMatch.size());

        auto params = this->getEmbeddingParameters();
        auto confs = RDKit::DGeomHelpers::EmbedMultipleConfs(*mol, numConfs, params);

        std::vector<unsigned> confIds;
        for (auto const confId : confs) {
            confIds.emplace_back(confId);
        }

        std::vector<std::pair<int, double>> result;
        RDKit::MMFF::MMFFOptimizeMoleculeConfs(*mol, result, m_threads);

        for (auto const confId : confs) {
            auto score = RDKit::MolAlign::alignMol(*mol, *m_core.ref, confId, 0, &molMatch);
            spdlog::debug("aligned conformer {} with score {}", confId, score);
        }

        assert(mol->getNumConformers() == numConfs);
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    RDKit::SubstructMatchParameters ConformerEmbedder::getSubstructMatchParameters() const {
        RDKit::SubstructMatchParameters params;
        params.uniquify = false;
        params.useChirality = false;
        params.useQueryQueryMatches = false;
        params.maxMatches = 1000;
        params.numThreads = m_threads;
        return params;
    }

    RDKit::DGeomHelpers::EmbedParameters ConformerEmbedder::getEmbeddingParameters() const {
        auto params = RDKit::DGeomHelpers::srETKDGv3;
        params.randomSeed = seed;
        params.numThreads = m_threads;
        params.optimizerForceTol = forceTol;
        params.clearConfs = false;
        return params;
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    std::vector<multialign::PoseID> ConformerEmbedder::generateNewPosesForAssemblyLigand(
        const multialign::Ligand &worstLigand, const multialign::LigandVector &targets,
        const std::unordered_map<multialign::LigandID, multialign::PoseID> &conformerIDs,
        const core::PairwiseMCSMap &pairwiseStrictMCSMap, const core::PairwiseMCSMap &pairwiseRelaxedMCSMap) {
        // TODO maybe sanitize mols?
        std::vector<unsigned> newIds;
        RDKit::ROMol *ligandMol = (RDKit::ROMol *)worstLigand.getMoleculePtr();

        for (const multialign::Ligand &target : targets) {
            // find mcs
            const multialign::LigandID targetID = target.getID();
            if (conformerIDs.count(targetID) == 0) {
                continue;
            }

            const multialign::PoseID targetConformerID = conformerIDs.at(targetID);
            const RDKit::ROMol targetMol = target.getMolecule();
            RDKit::Conformer targetConformer;
            try {
                targetConformer = targetMol.getConformer(targetConformerID);
            } catch (std::runtime_error e) {
                spdlog::error(e.what());
            }

            RDKit::MatchVectType ligandMatchRelaxed, targetMatchRelaxed, ligandMatchStrict, targetMatchStrict;
            // RDKit::MatchVectType smallerIDLigandMatch, largerIDLigandMatch;
            std::string mcsStringRelaxed, mcsStringStrict;
            multialign::LigandPair ligandPair {worstLigand.getID(), targetID};

            // since mcs maps are accessed via ligand pair, i.e. smaller id first, we have to check in which
            // order ligand and target are.
            if (worstLigand.getID() < targetID) {
                std::tie(ligandMatchRelaxed, targetMatchRelaxed, mcsStringRelaxed)
                    = pairwiseRelaxedMCSMap.at(ligandPair);
                std::tie(ligandMatchStrict, targetMatchStrict, mcsStringStrict) = pairwiseStrictMCSMap.at(ligandPair);
            } else {
                std::tie(targetMatchRelaxed, ligandMatchRelaxed, mcsStringRelaxed)
                    = pairwiseRelaxedMCSMap.at(ligandPair);
                std::tie(targetMatchStrict, ligandMatchStrict, mcsStringStrict) = pairwiseStrictMCSMap.at(ligandPair);
            }

            CoreAtomMapping ligandMcsCoords;
            RDKit::DGeomHelpers::EmbedParameters params = get_embed_params_for_optimizer_generation();
            int addedID = -1;

            // try relaxed mcs first
            if (!ligandMatchRelaxed.empty() && !targetMatchRelaxed.empty()) {
                spdlog::debug("trying relaxed substructure approach.");
                ligandMcsCoords = getLigandMcsAtomCoordsFromTargetMatch(targetConformer.getPositions(),
                                                                        ligandMatchRelaxed, targetMatchRelaxed);
                params.coordMap = &ligandMcsCoords;
                try {
                    addedID = RDKit::DGeomHelpers::EmbedMolecule(*ligandMol, params);
                } catch (const std::runtime_error &e) {
                    spdlog::debug(e.what());
                }
            }

            // if relaxed mcs params didnt yield valid embedding, reattempt with strict mcs.
            if (addedID < 0 && !ligandMatchStrict.empty() && !targetMatchStrict.empty()) {
                spdlog::debug("flexible approach failed. Trying strict approach.");
                ligandMcsCoords = getLigandMcsAtomCoordsFromTargetMatch(targetConformer.getPositions(),
                                                                        ligandMatchStrict, targetMatchStrict);
                params.coordMap = &ligandMcsCoords;
                try {
                    addedID = RDKit::DGeomHelpers::EmbedMolecule(*ligandMol, params);
                } catch (const std::runtime_error &e) {
                    spdlog::debug(e.what());
                }
            }
            if (addedID < 0) {
                spdlog::debug("strict mcs confgen failed. mcs: {}, target: {}, mol {}", mcsStringStrict,
                              RDKit::MolToSmiles(*worstLigand.getMoleculePtr()), RDKit::MolToSmiles(targetMol));
                spdlog::debug("target conformer {}/{}: no viable pose generated.", targetID, targets.size());
                continue;
            }
            spdlog::debug("target conformer {}/{}: generated valid pose.", targetID, targets.size());
            const unsigned addedIDUnsigned = static_cast<unsigned>(addedID);
            newIds.push_back(addedIDUnsigned);
        }
        return newIds;
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    CoreAtomMapping ConformerEmbedder::getLigandMcsAtomCoordsFromTargetMatch(
        const RDGeom::POINT3D_VECT &targetCoords, const RDKit::MatchVectType &ligandMcsMatch,
        const RDKit::MatchVectType &targetMcsMatch) {
        CoreAtomMapping mcsCoords;
        for (const auto &[mcsAtomID, targetAtomID] : targetMcsMatch) {
            mcsCoords[mcsAtomID] = targetCoords.at(targetAtomID);
        }
        CoreAtomMapping ligandCoords;
        for (const auto &[mcsAtomID, ligandAtomID] : ligandMcsMatch) {
            ligandCoords[ligandAtomID] = mcsCoords[mcsAtomID];
        }
        return ligandCoords;
    }

}  // namespace coaler::embedder
