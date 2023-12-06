add_executable (coaler main.cpp)
add_subdirectory(lib)

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

install(TARGETS coaler
        RUNTIME DESTINATION ${CMAKE_BINARY_DIR}
        LIBRARY DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
