add_executable (aligner
        main.cpp
        SingleAligner.cpp
        SingleAligner.hpp
)


target_include_directories(aligner PUBLIC ${Boost_INCLUDE_DIRS})

target_link_libraries(aligner PUBLIC rdkit::RDKitSmilesParse)
target_link_libraries(aligner PUBLIC rdkit::RDKitGraphMol)
target_link_libraries(aligner PUBLIC rdkit::RDKitRDGeneral)
target_link_libraries(aligner PUBLIC rdkit::RDKitRDGeometryLib)
target_link_libraries(aligner PUBLIC rdkit::RDKitDataStructs)
target_link_libraries(aligner PUBLIC rdkit::RDKitRingDecomposerLib)

install(TARGETS aligner DESTINATION "."
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
)