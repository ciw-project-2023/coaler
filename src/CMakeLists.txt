set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_executable (coaler main.cpp)
add_subdirectory(lib)

set_property(TARGET coaler PROPERTY CXX_STANDARD 20)

target_link_libraries(coaler PUBLIC "-Wl,--disable-new-dtags")

find_package(spdlog REQUIRED)
find_package(Boost REQUIRED COMPONENTS
        program_options
)

target_include_directories(coaler PUBLIC ${Boost_INCLUDE_DIRS})
target_include_directories(coaler PUBLIC ${spdlog_INCLUDE_DIRS})

target_include_directories(coaler PRIVATE lib)
target_link_directories(coaler PRIVATE lib)

target_link_libraries(coaler PUBLIC spdlog::spdlog)
target_link_libraries(coaler PUBLIC coaler_lib)

find_package(OpenMP REQUIRED)
target_compile_options(coaler PRIVATE ${OpenMP_CXX_FLAGS})
target_link_libraries(coaler PUBLIC OpenMP::OpenMP_CXX)

install(TARGETS coaler
        RUNTIME DESTINATION ${CMAKE_BINARY_DIR}
        LIBRARY DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
