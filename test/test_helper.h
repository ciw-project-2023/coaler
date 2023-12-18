//
// Created by niklas on 12/9/23.
//
#include "GraphMol/RWMol.h"
#include "GraphMol/SmilesParse/SmilesParse.h"

#ifndef COALER_TEST_HELPER_H
#    define COALER_TEST_HELPER_H

namespace {
    RDKit::RWMOL_SPTR MolFromSmiles(const std::string &smiles) {
        return boost::make_shared<RDKit::RWMol>(*RDKit::SmilesToMol(smiles));
    }

    RDKit::ROMOL_SPTR ROMolFromSmiles(const std::string &smiles) {
        return boost::make_shared<RDKit::ROMol>(*RDKit::SmilesToMol(smiles));
    }
}  // namespace

#endif  // COALER_TEST_HELPER_H
