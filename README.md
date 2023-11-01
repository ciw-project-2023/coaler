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
$ conan install .. -s build_type=Debug
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
Just make sure that you have run conan install

Now you should be able to see the RDKit headers you linked under `External Libraries -> Header Search Paths -> include`

Have fun coding!
