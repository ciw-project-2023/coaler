find_package (Catch2 REQUIRED)
include_directories (${TEST_SOURCE_DIR}/src)
add_executable (Test test_runner.cpp)

target_link_libraries(Test PRIVATE Catch2::Catch2)

include(CTest)
include(Catch)
#catch_discover_tests(Test)


add_test (NAME MyTest COMMAND Test)