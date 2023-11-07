add_executable (aligner
        main.cpp
        SingleAligner.cpp SingleAligner.hpp parsers/FileParserSmiles.h
        parsers/FileNotFoundException.h parsers/FileParserSmiles.cpp)

find_package(spdlog REQUIRED)
find_package(rdkit REQUIRED)
find_package(Boost REQUIRED COMPONENTS program_options)

target_link_libraries(aligner PUBLIC "-Wl,--disable-new-dtags")

target_include_directories(aligner PUBLIC ${Boost_INCLUDE_DIRS})
target_include_directories(aligner PUBLIC ${spdlog_INCLUDE_DIRS})

set(rdkitlibs "RDKitAbbreviations" "RDKitAlignment" "RDKitCatalogs" "RDKitChemicalFeatures" "RDKitChemReactions" "RDKitChemTransforms"
  "RDKitCIPLabeler" "RDKitcoordgen" "RDKitDataStructs" "RDKitDepictor" "RDKitDeprotect" "RDKitDescriptors" "RDKitDistGeometry" "RDKitDistGeomHelpers"
  "RDKitEigenSolvers" "RDKitFileParsers" "RDKitFilterCatalog" "RDKitFingerprints" "RDKitFMCS" "RDKitForceFieldHelpers" "RDKitForceField"
  "RDKitFragCatalog" "RDKitga" "RDKitGeneralizedSubstruct" "RDKitGenericGroups" "RDKitGraphMol" "RDKithc" "RDKitInfoTheory" "RDKitmaeparser"
  "RDKitMarvinParser" "RDKitMMPA" "RDKitMolAlign" "RDKitMolCatalog" "RDKitMolChemicalFeatures" "RDKitMolDraw2D" "RDKitMolEnumerator" "RDKitMolHash"
  "RDKitMolInterchange" "RDKitMolStandardize" "RDKitMolTransforms" "RDKitO3AAlign" "RDKitOptimizer" "RDKitPartialCharges" "RDKitRascalMCES"
  "RDKitRDGeneral" "RDKitRDGeometryLib" "RDKitRDStreams" "RDKitReducedGraphs" "RDKitRGroupDecomposition" "RDKitRingDecomposerLib" "RDKitScaffoldNetwork"
  "RDKitShapeHelpers" "RDKitSimDivPickers" "RDKitSLNParse" "RDKitSmilesParse" "RDKitSubgraphs" "RDKitSubstructLibrary" "RDKitSubstructMatch"
  "RDKitTautomerQuery" "RDKitTrajectory")

foreach(lib IN LISTS rdkitlibs)
    target_link_libraries(aligner PUBLIC rdkit::${lib})
endforeach()

target_link_libraries(aligner PUBLIC spdlog::spdlog)

install(TARGETS aligner DESTINATION "."
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)
