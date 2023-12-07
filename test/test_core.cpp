#include <GraphMol/RDKitBase.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <catch2/catch.hpp>
#include <string>

#include "coaler/core/Core.hpp"

using namespace coaler::core;

TEST_CASE("Core_constructor", "[core]") {
    RDKit::ROMol mol1 = *RDKit::SmilesToMol("c1cc(ccc1N2CCOCC2=O)N3C[C@@H](OC3=O)CNC(=O)c4ccc(s4)Cl");
    RDKit::ROMol mol2 = *RDKit::SmilesToMol("CC1(C(=O)N(C(=S)N1c2ccc(c(c2)F)C(=O)NC)c3ccc(c(c3)C(F)(F)F)C#N)C");
    std::vector<RDKit::RWMol> molVec;
    molVec.emplace_back(mol1);
    molVec.emplace_back(mol2);

    Core coreMCS(molVec, CoreType::MCS);
    RDKit::ROMol coreScaffoldMCS = coreMCS.getCore();
    CHECK(RDKit::MolToSmarts(coreScaffoldMCS) == "[#6]1:[#6]:[#6](:[#6]:[#6]:[#6]:1)-[#7](-[#6]-[#6]-[#6])-[#6]");

    Core coreMurcko(molVec, CoreType::Murcko);
    RDKit::ROMol coreScaffoldMurcko = coreMurcko.getCore();
    CHECK(RDKit::MolToSmarts(coreScaffoldMurcko) == "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1");
}