/*
 * Copyright 2023 CoAler Group, all rights reserved.
 */

#include "ConformerEmbedder.hpp"

#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/FMCS/FMCS.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/ForceFieldHelpers/UFF/UFF.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <spdlog/spdlog.h>

#include <utility>

#include "../core/Matcher.hpp"
#include "../multialign/Forward.hpp"
#include "../multialign/models/Forward.hpp"
#include "Forward.hpp"
#include "GraphMol/FMCS/FMCS.h"

const unsigned seed = 42;
const float forceTol = 0.0135;

namespace {
    auto get_mcs(const RDKit::MOL_SPTR_VECT &mols) {
        RDKit::MCSParameters mcsParams;
        RDKit::MCSAtomCompareParameters atomCompParams;
        atomCompParams.MatchChiralTag = true;
        atomCompParams.MatchFormalCharge = true;
        atomCompParams.MatchIsotope = false;
        atomCompParams.MatchValences = true;
        atomCompParams.RingMatchesRingOnly = true;
        atomCompParams.CompleteRingsOnly = true;
        mcsParams.AtomCompareParameters = atomCompParams;

        RDKit::MCSBondCompareParameters bondCompParams;
        bondCompParams.MatchStereo = true;
        bondCompParams.RingMatchesRingOnly = false;
        bondCompParams.CompleteRingsOnly = true;
        bondCompParams.MatchFusedRings = true;
        bondCompParams.MatchFusedRingsStrict = false;
        mcsParams.BondCompareParameters = bondCompParams;
        mcsParams.Timeout = 1;

        mcsParams.setMCSAtomTyperFromEnum(RDKit::AtomCompareElements);
        mcsParams.setMCSBondTyperFromEnum(RDKit::BondCompareOrderExact);

        RDKit::MCSResult const mcs = RDKit::findMCS(mols, &mcsParams);
        return mcs;
    }
}  // namespace

namespace coaler::embedder {
    ConformerEmbedder::ConformerEmbedder(const core::CoreResult &result, const int threads,
                                         const bool divideConformersByMatches)
        : m_core(result), m_threads(threads), m_divideConformersByMatches(divideConformersByMatches) {}

    void ConformerEmbedder::embedConformers(const RDKit::ROMOL_SPTR &mol, unsigned numConfs) {
        // firstMatch molecule and core
        RDKit::SubstructMatchParameters substructMatchParams;
        substructMatchParams.uniquify = false;
        substructMatchParams.useChirality = false;
        substructMatchParams.useQueryQueryMatches = false;
        substructMatchParams.maxMatches = 1000;
        substructMatchParams.numThreads = m_threads;

        auto matches = RDKit::SubstructMatch(*mol, *m_core.core, substructMatchParams);

        assert(!matches.empty());

        spdlog::debug("number of Core Matches: {}", matches.size());

        unsigned matchCounter = 0;
        for (auto const &match : matches) {
            auto params = this->getEmbeddingParameters();

            std::vector<int> confs;
            if (m_divideConformersByMatches) {
                auto numConfsMatch = (matchCounter < numConfs % matches.size()) ? (int(numConfs / matches.size()) + 1)
                                                                                : (int(numConfs / matches.size()));
                confs = RDKit::DGeomHelpers::EmbedMultipleConfs(*mol, numConfsMatch, params);

                matchCounter++;
            } else {
                confs = RDKit::DGeomHelpers::EmbedMultipleConfs(*mol, numConfs, params);
            }

            std::vector<unsigned> confIds;
            for (auto const confId : confs) {
                confIds.emplace_back(confId);
            }

            std::vector<std::pair<int, double>> result;
            RDKit::UFF::UFFOptimizeMoleculeConfs(*mol, result, m_threads);

            // This is somewhat unintuitive:
            // We use a reference molecule (m_core.ref) to have a "real" conformer, as we cannot generate
            // chemically sensible conformers from m_core.core, as it contains smarts queries. So we get the
            // match of the query to the reference in m_core.coreToRef {'id_in_core': 'id_in_ref'}. Now we have to make
            // a list for the alignment of [(id_in_mol, id_in_ref)] because the alignment wants the atom mapping in the
            // opposite order than we get from the substruct matching (i.e. we get (queryId, molId) from substruct
            // and have to provide (molId, queryId) to the alignment.
            RDKit::MatchVectType matchReverse;
            for (const auto &[queryId, molId] : match) {
                matchReverse.emplace_back(std::make_pair(molId, m_core.coreToRef.at(queryId)));
            }

            for (auto const confId : confs) {
                auto score = RDKit::MolAlign::alignMol(*mol, *m_core.ref, confId, 0, &matchReverse);
                spdlog::debug("aligned conformer {} with score {}", confId, score);
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
        params.randomSeed = seed;
        params.numThreads = m_threads;
        params.optimizerForceTol = forceTol;
        params.clearConfs = false;
        return params;
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    std::vector<multialign::PoseID> ConformerEmbedder::generateNewPosesForAssemblyLigand(
        RDKit::ROMol *worstLigandMol, const multialign::LigandVector &targets,
        const std::unordered_map<multialign::LigandID, multialign::PoseID> &conformerIDs) {
        std::vector<unsigned> newIds;
        for (const multialign::Ligand &target : targets) {
            // find mcs
            const multialign::LigandID targetID = target.getID();
            const multialign::PoseID targetConformerID = conformerIDs.at(targetID);
            const RDKit::ROMol targetMol = target.getMolecule();
            const RDKit::Conformer &targetConformer = targetMol.getConformer(targetConformerID);

            std::vector<RDKit::ROMOL_SPTR> mols = {
                boost::make_shared<RDKit::ROMol>(*worstLigandMol),
                boost::make_shared<RDKit::ROMol>(target.getMolecule()),
            };

            // const core::CoreResult mcsResult = matcher.coaler::core::Matcher::calculateCoreMcs(mols).value();

            auto mcsResult = get_mcs(mols);

            RDKit::SubstructMatchParameters substructMatchParams;
            substructMatchParams.uniquify = true;
            substructMatchParams.useChirality = true;
            substructMatchParams.useQueryQueryMatches = false;
            substructMatchParams.maxMatches = 1;
            substructMatchParams.numThreads = 1;

            std::vector<RDKit::MatchVectType> ligandMatches;
            RDKit::SubstructMatch(*worstLigandMol, *mcsResult.QueryMol, ligandMatches, substructMatchParams.uniquify, false,
                                  substructMatchParams.useChirality, substructMatchParams.useQueryQueryMatches,
                                  substructMatchParams.maxMatches, substructMatchParams.numThreads);

            std::vector<RDKit::MatchVectType> targetMatches;
            RDKit::SubstructMatch(targetMol, *mcsResult.QueryMol, targetMatches, substructMatchParams.uniquify, false,
                                  substructMatchParams.useChirality, substructMatchParams.useQueryQueryMatches,
                                  substructMatchParams.maxMatches, substructMatchParams.numThreads);

            assert(!targetMatches.empty() && !ligandMatches.empty());
            const RDKit::MatchVectType ligandMatch = ligandMatches.at(0);
            const RDKit::MatchVectType targetMatch = targetMatches.at(0);

            const CoreAtomMapping ligandMcsCoords
                = getLigandMcsAtomCoordsFromTargetMatch(targetConformer.getPositions(), ligandMatch, targetMatch);

            RDKit::DGeomHelpers::EmbedParameters params;
            params = RDKit::DGeomHelpers::srETKDGv3;
            params.optimizerForceTol = forceTol;
            params.randomSeed = seed;
            params.coordMap = &ligandMcsCoords;
            params.useRandomCoords = true;
            params.numThreads = 1;
            params.clearConfs = false;
            int addedID = 0;
            try {
                addedID = RDKit::DGeomHelpers::EmbedMolecule(*worstLigandMol, params);
            } catch (const std::runtime_error& e) {
                spdlog::warn(e.what());
                addedID = -1;
            }
            if (addedID < 0) {
                spdlog::warn("unable to generate pose for: ligand: {}, target: {}, mcs smarts: {}",
                             RDKit::MolToSmiles(*worstLigandMol),
                             RDKit::MolToSmiles(target.getMolecule()),
                             mcsResult.SmartsString);
                continue;
            }
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
