# TestRdkit

I am using Clion as IDE

## Install with Conan
Install conan: https://conan.io/downloads

### Clone rdkit 
https://github.com/ciw-project-2023/rdkit

Checkout `conan` Branch (`git checkout conan`)

### Run conan install + create

```bash
$ cd rdkit
$ conan create . -s build_type=Debug --build missing
```

This can take a loooong while (for very fast PCs under 10 minutes)

### Get your client repo
For example purpose:
https://github.com/ciw-project-2023/RDKitSkeletonCPP

### Add what you need from rdkit in your CMakeLists.txt
```sh
$ cd RDKitSkeletonCPP
$ conan install . -s build_type=Debug
$ $EDITOR CMakeLists.txt
```

```cmake
cmake_minimum_required(VERSION 3.15..3.27)
project(aligner CXX)

find_package(rdkit REQUIRED)

add_executable(aligner src/main.cpp src/SingleAligner.cpp src/SingleAligner.hpp)

target_link_libraries(aligner PUBLIC rdkit::RDKitSmilesParse)
target_link_libraries(aligner PUBLIC rdkit::RDKitGraphMol)
target_link_libraries(aligner PUBLIC rdkit::RDKitRDGeneral)
# add other libraries you need here

install(TARGETS aligner DESTINATION "."
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        )
```

The full list of linkable libs is:
-  RDKitAbbreviations
-  RDKitAlignment
-  RDKitCatalogs
-  RDKitChemicalFeatures
-  RDKitChemReactions
-  RDKitChemTransforms
-  RDKitCIPLabeler
-  RDKitcoordgen
-  RDKitDataStructs
-  RDKitDepictor
-  RDKitDeprotect
-  RDKitDescriptors
-  RDKitDistGeometry
-  RDKitDistGeomHelpers
-  RDKitEigenSolvers
-  RDKitFileParsers
-  RDKitFilterCatalog
-  RDKitFingerprints
-  RDKitFMCS
-  RDKitForceFieldHelpers
-  RDKitForceField
-  RDKitFragCatalog
-  RDKitga
-  RDKitGeneralizedSubstruct
-  RDKitGenericGroups
-  RDKitGraphMol
-  RDKithc
-  RDKitInfoTheory
-  RDKitmaeparser
-  RDKitMarvinParser
-  RDKitMMPA
-  RDKitMolAlign
-  RDKitMolCatalog
-  RDKitMolChemicalFeatures
-  RDKitMolDraw2D
-  RDKitMolEnumerator
-  RDKitMolHash
-  RDKitMolInterchange
-  RDKitMolStandardize
-  RDKitMolTransforms
-  RDKitO3AAlign
-  RDKitOptimizer
-  RDKitPartialCharges
-  RDKitRascalMCES
-  RDKitRDGeneral
-  RDKitRDGeometryLib
-  RDKitRDStreams
-  RDKitReducedGraphs
-  RDKitRGroupDecomposition
-  RDKitRingDecomposerLib
-  RDKitScaffoldNetwork
-  RDKitShapeHelpers
-  RDKitSimDivPickers
-  RDKitSLNParse
-  RDKitSmilesParse
-  RDKitSubgraphs
-  RDKitSubstructLibrary
-  RDKitSubstructMatch
-  RDKitTautomerQuery
-  RDKitTrajectory

### Have fun
In your `*.cpp` or `*.h` or `*.hpp` you can now include your RDKit headers:
```c++
#include <cstdint>
#include <iostream>

#include <GraphMol/SmilesParse/SmilesParse.h>

#include "SingleAligner.hpp"

int main(int argc, char *argv[]) {
RDKit::RWMol* mol_a = RDKit::SmilesToMol("CCCO");
RDKit::RWMol* mol_b = RDKit::SmilesToMol("CCCN");

ciw::SingleAligner single_aligner;

//single_aligner.algin_molecules(mol_a->, mol_b);

std::cout << "RDKit is linked :)" << std::endl;
}
```

### Run build
run
```bash
$ cd RDKitSkeletonCPP
$ conan build . -s build_type=Debug
```

Here you should get an output containing lines like this:
```
======== Calling build() ========
conanfile.py (aligner/0.0.1): Calling build()
conanfile.py (aligner/0.0.1): Running CMake.configure()
conanfile.py (aligner/0.0.1): RUN: cmake -G "Unix Makefiles" -DCMAKE_TOOLCHAIN_FILE="/home/niklas/projects/uni/RDKitSkeletonCPP/build/Debug/generators/conan_toolchain.cmake" -DCMAKE_INSTALL_PREFIX="/home/niklas/projects/uni/RDKitSkeletonCPP" -DCMAKE_POLICY_DEFAULT_CMP0091="NEW" -DCMAKE_BUILD_TYPE="Debug" "/home/niklas/projects/uni/RDKitSkeletonCPP"
```

Copy this part (it will be different for on your PC): 
```
-DCMAKE_TOOLCHAIN_FILE="/home/niklas/projects/uni/RDKitSkeletonCPP/build/Debug/generators/conan_toolchain.cmake" -DCMAKE_INSTALL_PREFIX="/home/niklas/projects/uni/RDKitSkeletonCPP" `
```

### Setup Clion
1. Open your Project in Clion and tell it that it is a CMake project
2. In Clion Goto Settings->Build,Execution...->CMAKE and add the thing you just copied to `CMake Options`.
3. You'll get a CMake plugin tab in your plugins -> Run "Reload CMake"
4. If your build worked, this should work without error. 

Now you should be able to see the RDKit headers you linked under `External Libraries -> Header Search Paths -> include`

Have fun coding!