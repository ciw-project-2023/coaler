cmake_minimum_required(VERSION 3.26)
project(TestRDkit)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


find_package(RDKit REQUIRED)
find_package(Boost COMPONENTS timer system REQUIRED)

add_executable(TestRDkit src/main.cpp)

target_link_libraries(TestRDkit RDKit::SmilesParse RDKit::MolHash
        Boost::timer)
