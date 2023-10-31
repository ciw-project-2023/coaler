add_executable (aligner
        main.cpp
        SingleAligner.cpp
)

target_link_libraries(aligner PUBLIC rdkit::RDKitSmilesParse)
target_link_libraries(aligner PUBLIC rdkit::RDKitGraphMol)
target_link_libraries(aligner PUBLIC rdkit::RDKitRDGeneral)

install(TARGETS aligner DESTINATION "."
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
)