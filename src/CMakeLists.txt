add_executable (aligner
        main.cpp
        SingleAligner.cpp
        SingleAligner.hpp
        multialign/PairwiseAlignment.cpp
        multialign/PairwiseAlignment.hpp
        multialign/PoseRegister.cpp
        multialign/PoseRegister.hpp
        multialign/Forward.hpp
        multialign/Ligand.cpp
        multialign/Ligand.hpp
        multialign/PosePair.cpp
        multialign/PosePair.hpp multialign/MultiAligner.cpp multialign/MultiAligner.hpp multialign/LigandAlignmentAssembly.cpp multialign/LigandAlignmentAssembly.hpp multialign/PoseRegisterBuilder.cpp multialign/PoseRegisterBuilder.hpp multialign/LigandPair.cpp multialign/LigandPair.hpp
)

target_link_libraries(aligner PUBLIC "-Wl,--disable-new-dtags")

target_include_directories(aligner PUBLIC ${Boost_INCLUDE_DIRS})

set(rdkitlibs "RDKitAbbreviations" "RDKitAlignment" "RDKitCatalogs" "RDKitChemicalFeatures" "RDKitChemReactions" "RDKitChemTransforms"
  "RDKitCIPLabeler" "RDKitcoordgen" "RDKitDataStructs" "RDKitDepictor" "RDKitDeprotect" "RDKitDescriptors" "RDKitDistGeometry" "RDKitDistGeomHelpers"
  "RDKitEigenSolvers" "RDKitFileParsers" "RDKitFilterCatalog" "RDKitFingerprints" "RDKitFMCS" "RDKitForceFieldHelpers" "RDKitForceField"
  "RDKitFragCatalog" "RDKitga" "RDKitGeneralizedSubstruct" "RDKitGenericGroups" "RDKitGraphMol" "RDKithc" "RDKitInfoTheory" "RDKitmaeparser"
  "RDKitMarvinParser" "RDKitMMPA" "RDKitMolAlign" "RDKitMolCatalog" "RDKitMolChemicalFeatures" "RDKitMolDraw2D" "RDKitMolEnumerator" "RDKitMolHash"
  "RDKitMolInterchange" "RDKitMolStandardize" "RDKitMolTransforms" "RDKitO3AAlign" "RDKitOptimizer" "RDKitPartialCharges" "RDKitRascalMCES"
  "RDKitRDGeneral" "RDKitRDGeometryLib" "RDKitRDStreams" "RDKitReducedGraphs" "RDKitRGroupDecomposition" "RDKitRingDecomposerLib" "RDKitScaffoldNetwork"
  "RDKitShapeHelpers" "RDKitSimDivPickers" "RDKitSLNParse" "RDKitSmilesParse" "RDKitSubgraphs" "RDKitSubstructLibrary" "RDKitSubstructMatch"
  "RDKitTautomerQuery" "RDKitTrajectory")

foreach(lib IN LISTS rdkitlibs)
    target_link_libraries(aligner PUBLIC rdkit::${lib})
endforeach()


install(TARGETS ${rdkit_LIBS} RUNTIME DESTINATION bin)

install(TARGETS aligner DESTINATION "."
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
)