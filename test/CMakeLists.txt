find_package(Catch2 REQUIRED)
find_package(Boost REQUIRED)

include_directories (${TEST_SOURCE_DIR}/src)

file(GLOB TEST_SOURCES "*.hpp" "*.cpp")

add_executable (Test test_runner.cpp
        pose_register_tester.cpp
        ${CMAKE_SOURCE_DIR}/src/multialign/BasicClasses/PosePair.cpp
        ${CMAKE_SOURCE_DIR}/src/multialign/PoseRegister.cpp
        ${CMAKE_SOURCE_DIR}/src/multialign/BasicClasses/Ligand.cpp
        ${CMAKE_SOURCE_DIR}/src/multialign/BasicClasses/LigandPair.cpp
        ligand_tester.cpp multi_aligner_tester.cpp
        ligand_pair_tester.cpp
        pose_pair_tester.cpp
        ligand_alignment_assembly_tester.cpp
        pose_register_builder_tester.cpp)


target_include_directories(Test PUBLIC ${Boost_INCLUDE_DIRS})

target_link_libraries(Test PUBLIC "-Wl,--disable-new-dtags")

target_link_libraries(Test PRIVATE Catch2::Catch2)
target_link_libraries(Test PUBLIC rdkit::RDKitSmilesParse)
target_link_libraries(Test PUBLIC rdkit::RDKitGraphMol)
target_link_libraries(Test PUBLIC rdkit::RDKitRDGeneral)
target_link_libraries(Test PUBLIC rdkit::RDKitRDGeometryLib)
target_link_libraries(Test PUBLIC rdkit::RDKitDataStructs)
target_link_libraries(Test PUBLIC rdkit::RDKitRingDecomposerLib)

include(CTest)
include(Catch)
#catch_discover_tests(Test)

add_test (NAME MyTest COMMAND Test)