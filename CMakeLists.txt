cmake_minimum_required(VERSION 3.15..3.27)
project(aligner CXX)

find_package(rdkit REQUIRED)

add_executable(aligner src/main.cpp src/SingleAligner.cpp src/SingleAligner.hpp)

target_link_libraries(aligner PUBLIC rdkit::RDKitSmilesParse)
target_link_libraries(aligner PUBLIC rdkit::RDKitGraphMol)
target_link_libraries(aligner PUBLIC rdkit::RDKitRDGeneral)

install(TARGETS aligner DESTINATION "."
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        )
