# TestRdkit

I am using Clion as IDE

## Install with Conan
Install conan: https://conan.io/downloads

### Fetch the rdkit package from Artifactory:

Credentials are:
```
username: ciw
password: VZKhzh2v5nCnijAS3A8R
```

```bash
$ conan remote add doc https://monkfish-app-qfnky.ondigitalocean.app
$ conan remote login doc
```

You can then proceed to the conan install step and it should pull rdkit from the server

### Alternatively: Run conan install + create

Clone https://github.com/ciw-project-2023/rdkit

Checkout `conan` Branch (`git checkout conan`)

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
$ conan install . -s build_type=Debug --build missing
```

### Run build
run
```bash
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
3. Change "Build directory" from `cmake-build-debug` to `build` (Conan uses this one)
4. You'll get a CMake plugin tab in your plugins -> Run "Reload CMake"
5. If your build worked, this should work without error. 

Now you should be able to see the RDKit headers you linked under `External Libraries -> Header Search Paths -> include`
