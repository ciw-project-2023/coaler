file(GLOB_RECURSE library_src CONFIGURE_DEPENDS *.cpp *.hpp *.h)

find_package(OpenMP REQUIRED)
find_package(rdkit REQUIRED)
find_package(Boost REQUIRED)
find_package(spdlog REQUIRED)

set(rdkitlibs "RDKitAbbreviations" "RDKitAlignment" "RDKitCatalogs" "RDKitChemicalFeatures" "RDKitChemReactions" "RDKitChemTransforms"
        "RDKitCIPLabeler" "RDKitcoordgen" "RDKitDataStructs" "RDKitDepictor" "RDKitDeprotect" "RDKitDescriptors" "RDKitDistGeometry" "RDKitDistGeomHelpers"
        "RDKitEigenSolvers" "RDKitFileParsers" "RDKitFilterCatalog" "RDKitFingerprints" "RDKitFMCS" "RDKitForceFieldHelpers" "RDKitForceField"
        "RDKitFragCatalog" "RDKitga" "RDKitGeneralizedSubstruct" "RDKitGenericGroups" "RDKitGraphMol" "RDKithc" "RDKitInfoTheory" "RDKitmaeparser"
        "RDKitMarvinParser" "RDKitMMPA" "RDKitMolAlign" "RDKitMolCatalog" "RDKitMolChemicalFeatures" "RDKitMolDraw2D" "RDKitMolEnumerator" "RDKitMolHash"
        "RDKitMolInterchange" "RDKitMolStandardize" "RDKitMolTransforms" "RDKitO3AAlign" "RDKitOptimizer" "RDKitPartialCharges" "RDKitRascalMCES"
        "RDKitRDGeneral" "RDKitRDGeometryLib" "RDKitRDStreams" "RDKitReducedGraphs" "RDKitRGroupDecomposition" "RDKitRingDecomposerLib" "RDKitScaffoldNetwork"
        "RDKitShapeHelpers" "RDKitSimDivPickers" "RDKitSLNParse" "RDKitSmilesParse" "RDKitSubgraphs" "RDKitSubstructLibrary" "RDKitSubstructMatch"
        "RDKitTautomerQuery" "RDKitTrajectory")


add_library(coaler_lib STATIC ${library_src} coaler/embedder/BucketDistributor.cpp coaler/embedder/BucketDistributor.hpp)

target_include_directories(coaler_lib PUBLIC ${rdkit_INCLUDE_DIRS})
target_include_directories(coaler_lib PUBLIC ${Boost_INCLUDE_DIRS})
target_include_directories(coaler_lib PUBLIC ${spdlog_INCLUDE_DIRS})

foreach(lib IN LISTS rdkitlibs)
    target_link_libraries(coaler_lib PUBLIC rdkit::${lib})
endforeach()

target_link_libraries(coaler_lib PUBLIC OpenMP::OpenMP_CXX)
target_link_libraries(coaler_lib PUBLIC spdlog::spdlog)
target_link_libraries(coaler_lib PUBLIC ${Boost_LIBRARIES})
