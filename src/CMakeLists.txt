add_executable (coaler
        main.cpp
        SingleAligner.cpp
        parsers/FileParserSmiles.cpp
)

find_package(spdlog REQUIRED)
find_package(rdkit REQUIRED)
find_package(Boost REQUIRED COMPONENTS program_options)

target_link_libraries(coaler PUBLIC "-Wl,--disable-new-dtags")

target_include_directories(coaler PUBLIC ${Boost_INCLUDE_DIRS})
target_include_directories(coaler PUBLIC ${spdlog_INCLUDE_DIRS})

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
    target_link_libraries(coaler PUBLIC rdkit::${lib})
endforeach()

target_link_libraries(coaler PUBLIC spdlog::spdlog)

install(TARGETS coaler DESTINATION "."
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)
