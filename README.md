# TestRdkit

I am using Clion as IDE

# Install
Install RDKit:

```commandline
conda create -c conda-forge -n my-rdkit-env rdkit
conda activate my-rdkit-env
```

Install Eigen3
```commandline
sudo apt install libeigen3-dev
```

# Error

error: ‘uint64_t’ in namespace ‘std’ does not name a type; did you mean ‘wint_t’?

Solution:
```c++
#include <cstdint> // before any RDKit includes
```

# CLion settings

First run:
```commandline
conda activate my-rdkit-env
cd PathToProject/RDKitSkeletonCPP
cmake .
```

## After that you can set the CLion settings:

Generator: Let CMake decide

CMake options:
```commandline
-DCMAKE_BUILD_TYPE=Debug
```

Now reload your CMakeLists.txt and the libraries should be found.
