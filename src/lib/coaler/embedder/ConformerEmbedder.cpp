#include "ConformerEmbedder.hpp"

#include <GraphMol/Atom.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/FMCS/FMCS.h>
#include <GraphMol/ForceFieldHelpers/UFF/UFF.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <spdlog/spdlog.h>

#include <utility>

const unsigned SEED = 42;
const float FORCE_TOL = 0.0135;
const unsigned BRUTEFORCE_CONFS = 500;

namespace {
    RDKit::SubstructMatchParameters get_optimizer_substruct_params() {
        RDKit::SubstructMatchParameters substructMatchParams;
        substructMatchParams.uniquify = true;
        substructMatchParams.useChirality = false;
        substructMatchParams.useQueryQueryMatches = false;
        substructMatchParams.maxMatches = 1;
        substructMatchParams.numThreads = 1;
        substructMatchParams.useEnhancedStereo = false;
        return substructMatchParams;
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    RDKit::DGeomHelpers::EmbedParameters get_embed_params_for_optimizer_generation() {
        RDKit::DGeomHelpers::EmbedParameters params;
        params = RDKit::DGeomHelpers::srETKDGv3;
        params.optimizerForceTol = FORCE_TOL;
        params.randomSeed = SEED;
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

    void ConformerEmbedder::embedConformers(const RDKit::ROMOL_SPTR &mol, unsigned numConfs) {
        spdlog::debug("embedding on: {}", RDKit::MolToSmiles(*mol));
        auto params = this->getEmbeddingParameters();

        if (RDKit::MolToSmiles(*mol) == RDKit::MolToSmiles(*m_core.ref)) {
            RDKit::DGeomHelpers::EmbedMultipleConfs(*mol, numConfs, params);

            return;
        }

        RDKit::MOL_SPTR_VECT pair = {mol, m_core.ref};

        auto mcsParams = core::Matcher::getRelaxedMCSParams();
        mcsParams.InitialSeed = RDKit::MolToSmarts(*m_core.core);

        auto mcs = RDKit::findMCS(pair, &mcsParams);
        if (mcs.QueryMol == nullptr) {
            spdlog::error("failed to find MCS for {} and {}", RDKit::MolToSmiles(*mol),
                          RDKit::MolToSmiles(*m_core.ref));
            return;
        }

        // firstMatch molecule and core
        RDKit::SubstructMatchParameters substructMatchParams;
        substructMatchParams.uniquify = false;
        substructMatchParams.useChirality = true;
        substructMatchParams.useQueryQueryMatches = true;
        substructMatchParams.useEnhancedStereo = true;
        substructMatchParams.maxMatches = 1000;
        substructMatchParams.numThreads = m_threads;

        assert(mcs.QueryMol.get()->getNumAtoms() >= m_core.core.get()->getNumAtoms());

        auto matches = RDKit::SubstructMatch(*mol, *mcs.QueryMol, substructMatchParams);
        spdlog::debug("number of mol matches: {}", matches.size());

        auto refMatch = RDKit::SubstructMatch(*m_core.ref, *mcs.QueryMol, substructMatchParams).at(0);
        std::unordered_map<int, int> refMatchMap;
        for (const auto &[queryId, molId] : refMatch) {
            refMatchMap[queryId] = molId;
        }

        unsigned matchCounter = 0;
        if (m_divideConformersByMatches && (numConfs / (float)matches.size()) < 5) {
            spdlog::warn("adding more conformers to get at least 5 per match");

            numConfs = matches.size() * 5;
        }

        for (auto const &match : matches) {

            std::vector<int> confs;
            if (m_divideConformersByMatches) {
                auto numConfsMatch = (matchCounter < numConfs % matches.size()) ? (int(numConfs / matches.size()) + 1)
                                                                                : (int(numConfs / matches.size()));
                confs = RDKit::DGeomHelpers::EmbedMultipleConfs(*mol, numConfsMatch, params);

                matchCounter++;
            } else {
                confs = RDKit::DGeomHelpers::EmbedMultipleConfs(*mol, numConfs, params);
            }

            std::vector<std::pair<int, double>> result;
            RDKit::UFF::UFFOptimizeMoleculeConfs(*mol, result, m_threads);

            RDKit::MatchVectType atomMap;
            for (const auto &[queryId, molId] : match) {
                    atomMap.emplace_back(molId, refMatchMap.at(queryId));
            }

            for (auto const confId : confs) {
                try {
                    auto score = RDKit::MolAlign::alignMol(*mol, *m_core.ref, confId, 0, &atomMap);
                    spdlog::debug("aligned conformer {} with score {}", confId, score);
                } catch (RDKit::MolAlign::MolAlignException &e) {
                    spdlog::debug("failed to align conformer {} with {}, {}: {}", confId, RDKit::MolToSmiles(*mol), RDKit::MolToSmiles(*m_core.ref), e.what());
                }
            }
        }

        if (m_divideConformersByMatches) {
            assert(mol->getNumConformers() == numConfs);
        } else {
            assert(mol->getNumConformers() == numConfs * matches.size());
        }
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    RDKit::DGeomHelpers::EmbedParameters ConformerEmbedder::getEmbeddingParameters() const {
        auto params = RDKit::DGeomHelpers::srETKDGv3;
        params.randomSeed = SEED;
        params.numThreads = m_threads;
        params.optimizerForceTol = FORCE_TOL;
        params.clearConfs = false;
        return params;
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    std::vector<multialign::PoseID> ConformerEmbedder::generateNewPosesForAssemblyLigand(
        const multialign::Ligand &worstLigand, const multialign::LigandVector &targets,
        const std::unordered_map<multialign::LigandID, multialign::PoseID> &conformerIDs,
        const core::PairwiseMCSMap &pairwiseStrictMCSMap, const core::PairwiseMCSMap &pairwiseRelaxedMCSMap,
        bool enforceGeneration) {
        std::vector<unsigned> newIds;
        auto *ligandMol = (RDKit::ROMol *)worstLigand.getMoleculePtr();

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
                targetConformer = targetMol.getConformer(static_cast<int>(targetConformerID));
            } catch (std::runtime_error &e) {
                spdlog::error(e.what());
            }

            // Clang-tidy readability!
            RDKit::MatchVectType ligandMatchRelaxed;
            RDKit::MatchVectType targetMatchRelaxed;
            RDKit::MatchVectType ligandMatchStrict;
            RDKit::MatchVectType targetMatchStrict;
            // RDKit::MatchVectType smallerIDLigandMatch, largerIDLigandMatch;
            std::string mcsStringRelaxed;
            std::string mcsStringStrict;
            const multialign::LigandPair ligandPair(worstLigand.getID(), targetID);

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

            double relaxedMcsSizeFactor = (double)ligandMatchRelaxed.size() / ligandMol->getNumAtoms();
            if (!(relaxedMcsSizeFactor > 0.2 || enforceGeneration)) {
                spdlog::debug("skipped due to small mcs, {} / {} = {}", ligandMatchRelaxed.size(),
                              ligandMol->getNumAtoms(), relaxedMcsSizeFactor);
                continue;
            }

            // check whether mapped chiral atoms match. Otherwise the embedder will take very long and then fail.
            bool invalidChiral = false;
            for (const auto &[targetMcsAtomID, targetAtomID] : targetMatchRelaxed) {
                auto targetChiralTag = targetMol.getAtomWithIdx(targetAtomID)->getChiralTag();
                if (targetChiralTag != RDKit::Atom::CHI_UNSPECIFIED) {
                    int matchingLigandAtomID = -1;
                    for (const auto &[ligandMcsAtomID, ligandAtomID] : ligandMatchRelaxed) {
                        if (ligandMcsAtomID == targetMcsAtomID) {
                            matchingLigandAtomID = ligandAtomID;
                            break;
                        }
                    }
                    assert(matchingLigandAtomID != -1);
                    auto ligandChiralTag = ligandMol->getAtomWithIdx(matchingLigandAtomID)->getChiralTag();
                    // todo sometimes chiral vs no chiral is not caught here, not sure why.
                    // embedding often fails in this case
                    // CHI_UNSPECIFIED (=0) sollte eigentlich bei anderen CHI tags != auswerten, idk was da falsch läuft
                    if (ligandChiralTag != targetChiralTag) {
                        spdlog::debug(
                            "chirality mismatch: \n"
                            "{}\n{}",
                            RDKit::MolToSmiles(*worstLigand.getMoleculePtr()), RDKit::MolToSmiles(targetMol));
                        invalidChiral = true;
                        break;
                    }
                }
            }
            if (invalidChiral) {
                continue;
            }

            // try relaxed mcs first
            if (!ligandMatchRelaxed.empty() && !targetMatchRelaxed.empty()) {
                ligandMcsCoords = getLigandMcsAtomCoordsFromTargetMatch(targetConformer.getPositions(),
                                                                        ligandMatchRelaxed, targetMatchRelaxed);
                params.coordMap = &ligandMcsCoords;
                try {
                    auto start = std::chrono::high_resolution_clock::now();
                    addedID = RDKit::DGeomHelpers::EmbedMolecule(*ligandMol, params);
                    auto end = std::chrono::high_resolution_clock::now();
                    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
                    if (duration > 5000) {
                        spdlog::debug("relaxed mcs confgen took {} ms", duration);
                        spdlog::debug("mol1: {} \nmol2: {}\n mcs: {}",
                                      RDKit::MolToSmiles(*worstLigand.getMoleculePtr()), RDKit::MolToSmiles(targetMol),
                                      mcsStringRelaxed);
                        spdlog::debug("success: {}", addedID > 0 ? "true\n" : "false\n");
                    }
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
                    auto start = std::chrono::high_resolution_clock::now();
                    addedID = RDKit::DGeomHelpers::EmbedMolecule(*ligandMol, params);
                    auto end = std::chrono::high_resolution_clock::now();
                    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
                    if (duration > 5000) {
                        spdlog::debug("strict mcs confgen took {} ms", duration);
                        spdlog::debug("mol1: {} \nmol2: {}\n mcs: {}\n",
                                      RDKit::MolToSmiles(*worstLigand.getMoleculePtr()), RDKit::MolToSmiles(targetMol),
                                      mcsStringStrict);
                    }
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
            const auto addedIDUnsigned = static_cast<unsigned>(addedID);
            newIds.push_back(addedIDUnsigned);
        }
        return newIds;
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    std::vector<multialign::PoseID> ConformerEmbedder::generateNewPosesForAssemblyLigand(
        const multialign::Ligand &worstLigand) {
        std::vector<unsigned> newIds;
        std::vector<int> newIntIds;
        auto *ligandMol = (RDKit::ROMol *)worstLigand.getMoleculePtr();
        RDKit::DGeomHelpers::EmbedParameters params = get_embed_params_for_optimizer_generation();

        RDKit::MatchVectType ligandMatch;
        RDKit::MatchVectType targetMatch;

        // get atom coords for core structure
        const CoreAtomMapping coreCoords = getLigandMcsAtomCoordsFromTargetMatch(
            m_core.ref->getConformer(0).getPositions(), ligandMatch, targetMatch);
        params.coordMap = &coreCoords;
        params.numThreads = m_threads;

        // embed BRUTEFORCE_CONFS new conformers into ligand
        try {
            newIntIds = RDKit::DGeomHelpers::EmbedMultipleConfs(*ligandMol, BRUTEFORCE_CONFS, params);
        } catch (const std::runtime_error &e) {
            spdlog::debug(e.what());
        }

        // type casting
        for (auto newIntId : newIntIds) {
            unsigned newUnsignedId = static_cast<unsigned>(newIntId);
            newIds.push_back(newUnsignedId);
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
