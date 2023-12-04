find_package(Catch2 REQUIRED)
find_package(spdlog REQUIRED)
find_package(Boost REQUIRED)

file(GLOB test_src CONFIGURE_DEPENDS "*.hpp" "*.cpp")



add_executable (Test ${test_src} conformer_generator_tester.cpp)

target_link_libraries(Test PUBLIC "-Wl,--disable-new-dtags")

# Include additional headers here
target_include_directories(Test PUBLIC
        ${Boost_INCLUDE_DIRS}
        ${CMAKE_SOURCE_DIR}/src/lib/
)

target_link_libraries(Test PRIVATE Catch2::Catch2)
target_link_libraries(Test PUBLIC spdlog::spdlog)
target_link_libraries(Test PUBLIC coaler_lib)

set(TEST_DATA_FILES AID_1806504.smi)
add_test(NAME MyTest COMMAND Test)
