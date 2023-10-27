cmake_minimum_required(VERSION 3.26)
project(TestRDkit)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


find_package(RDKit REQUIRED)
find_package(Boost COMPONENTS timer system REQUIRED)

add_executable(TestRDkit src/main.cpp)

target_link_libraries(TestRDkit
        RDKit::SmilesParse
        RDKit::MolHash
        RDKit::MolAlign
        RDKit::ga
        RDKit::hc
        RDKit::FMCS
        RDKit::MMPA
        RDKit::RDBoost
        RDKit::Catalogs
        RDKit::Depictor
        RDKit::GraphMol
        RDKit::O3AAlign
        RDKit::SLNParse
        RDKit::Alignment
        RDKit::Deprotect
        #RDKit::MolDraw2D Needs Cairo?
        RDKit::Optimizer
        RDKit::RDGeneral
        #RDKit::RDStreams needs Boost::iostreams
        RDKit::Subgraphs
        RDKit::CIPLabeler
        RDKit::ForceField
        RDKit::InfoTheory
        RDKit::MolCatalog
        RDKit::Trajectory
        RDKit::DataStructs
        RDKit::Descriptors
        RDKit::FileParsers
        RDKit::FragCatalog
        RDKit::SmilesParse
        RDKit::DistGeometry
        RDKit::EigenSolvers
        RDKit::Fingerprints
        RDKit::ShapeHelpers
        RDKit::Abbreviations
        RDKit::ChemReactions
        RDKit::FilterCatalog
        RDKit::GenericGroups
        RDKit::MolEnumerator
        RDKit::MolTransforms
        RDKit::RDGeometryLib
        RDKit::ReducedGraphs
        RDKit::SimDivPickers
        RDKit::TautomerQuery
        RDKit::ChemTransforms
        RDKit::MolInterchange
        RDKit::MolStandardize
        RDKit::PartialCharges
        RDKit::SubstructMatch
        RDKit::DistGeomHelpers
        RDKit::ScaffoldNetwork
        RDKit::ChemicalFeatures
        RDKit::SubstructLibrary
        RDKit::ForceFieldHelpers
        RDKit::MolHash
        Boost::timer)
