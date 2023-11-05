cmake_minimum_required(VERSION 3.15..3.27)
project(aligner CXX)

find_package(rdkit REQUIRED)

add_definitions(-DFMT_HEADER_ONLY)

add_subdirectory(src)
enable_testing()
add_subdirectory (test)
