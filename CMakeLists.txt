cmake_minimum_required(VERSION 3.15..3.27)
project(aligner CXX)

find_package(mytestlib REQUIRED)
find_package(rdkit REQUIRED)

add_executable(aligner src/aligner.cpp src/main.cpp)

target_link_libraries(aligner PUBLIC rdkit::RDKitGraphMol)
target_link_libraries(aligner PUBLIC rdkit::RDKitRDGeneral)
target_link_libraries(aligner PUBLIC mytestlib::mytest)
target_link_libraries(aligner PUBLIC mytestlib::second)

install(TARGETS aligner DESTINATION "."
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        )
