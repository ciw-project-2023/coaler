cmake_minimum_required(VERSION 3.26)
project(TestRDkit)

set(CMAKE_CXX_STANDARD 17)

find_package(RDKit REQUIRED)

add_executable(TestRDkit main.cpp)

target_link_libraries(TestRDkit RDKit::SmilesParse RDKit::MolHash)