find_package(Catch2 REQUIRED)
find_package(Boost REQUIRED)

include_directories (${TEST_SOURCE_DIR}/src)
add_executable (Test test_runner.cpp)

target_include_directories(Test PUBLIC ${Boost_INCLUDE_DIRS})

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