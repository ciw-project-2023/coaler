find_package(Catch2 REQUIRED)
find_package(spdlog REQUIRED)
find_package(Boost REQUIRED)

add_executable (Test
        test_runner.cpp
        test_single_aligner.cpp
        test_file_parser.cpp
        ${CMAKE_SOURCE_DIR}/src/singlealign/SingleAligner.cpp # TODO: restructure this
        ${CMAKE_SOURCE_DIR}/src/parser/FileParser.cpp
        pose_register_tester.cpp
        ${CMAKE_SOURCE_DIR}/src/multialign/BasicClasses/PosePair.cpp
        ${CMAKE_SOURCE_DIR}/src/multialign/PoseRegister.cpp
        ${CMAKE_SOURCE_DIR}/src/multialign/BasicClasses/Ligand.cpp
        ${CMAKE_SOURCE_DIR}/src/multialign/BasicClasses/LigandPair.cpp
        ${CMAKE_SOURCE_DIR}/src/multialign/PoseRegisterBuilder.cpp
        ligand_tester.cpp multi_aligner_tester.cpp
        ligand_pair_tester.cpp
        pose_pair_tester.cpp
        ligand_alignment_assembly_tester.cpp
        pose_register_builder_tester.cpp)

include_directories (${TEST_SOURCE_DIR}/src)


target_include_directories(Test PUBLIC ${Boost_INCLUDE_DIRS})

# Include additional headers here
target_include_directories(Test PUBLIC
        ${Boost_INCLUDE_DIRS}
        ${CMAKE_SOURCE_DIR}/src/
)

# All the dependencies get linked here
target_link_libraries(Test PUBLIC "-Wl,--disable-new-dtags")
target_link_libraries(Test PRIVATE Catch2::Catch2)
target_link_libraries(Test PUBLIC spdlog::spdlog)

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
    target_link_libraries(Test PUBLIC rdkit::${lib})
endforeach()

set(TEST_DATA_FILES AID_1806504.smi)
add_test(NAME MyTest COMMAND Test)
