# This is a skeleton project that uses RDKit

# Install
Install RDKit:

```commandline
sudo apt-get install python3-rdkit librdkit1 rdkit-data
```

Install Eigen3
```commandline
sudo apt install libeigen3-dev
```

# Error

I got this Error `error: ‘uint64_t’ in namespace ‘std’ does not name a type; did you mean ‘wint_t’?`

The solution was:
```c++
#include <cstdint> // before any RDKit includes
```

# CLion settings

Generator: Set to: `Let CMake decide

CMake options:
```commandline
-DCMAKE_BUILD_TYPE=Debug
```