# coaler

<img src="logo/coaler-removebg.png" width="150">

Core Aligner

### UPDATE YOUR SYSTEM

We have prebuild rdkit for the following setups:

```
[settings]
arch=x86_64
build_type=Debug
compiler=gcc
compiler.cppstd=gnu17
compiler.libcxx=libstdc++11
compiler.version=13
os=Linux

-- GlibC 2.36/2.38
```

And

```
[settings]
arch=x86_64
build_type=Debug
compiler=gcc
compiler.cppstd=gnu17
compiler.libcxx=libstdc++11
compiler.version=11
os=Linux

-- GlibC 2.34/2.36 (should work for both)
```

Please check your version you have installed locally and update them if you need to:

```bash
# Check g++ (compile.version)
g++ -v

# Check libc version
ldd --version

# Check cmake version 
cmake --version
```

### Install conan

<https://conan.io/downloads>

### Add your conan profile

Make conan detect the default profile:

```
conan profile detect
```

This will create `~/.conan2/profiles/default`.

Then for local development, you want a debug profile. So add change this in your
conan default profile.

~/.conan2/profiles/default

```
- build_type=Release
+ build_type=Debug
```

### Fetch the rdkit package from Artifactory

Credentials are:

```
username: ciw
password: VZKhzh2v5nCnijAS3A8R
```

```bash
conan remote add coaler http://server.conan.corealigner.de
```

```
conan remote login coaler
```

You can then proceed to the conan install step and it should pull rdkit from the server

### Clone the client Repo

<https://github.com/ciw-project-2023/coaler>

```sh
git clone https://github.com/ciw-project-2023/coaler && cd coaler
```

### Build Project

```bash
conan build . --build missing --update
```

To run the programm from the CMD:

```
./build/Debug/src/aligner
```

To run the tests from the CMD:

```
./build/Debug/test/Test
```

### Setup Clion

In the output of `conan build` you should get an output containing lines like this:

```
======== Calling build() ========
conanfile.py (aligner/0.0.1): Calling build()
conanfile.py (aligner/0.0.1): Running CMake.configure()
conanfile.py (aligner/0.0.1): RUN: cmake -G "Unix Makefiles" -DCMAKE_TOOLCHAIN_FILE="/home/niklas/projects/uni/coaler/build/Debug/generators/conan_toolchain.cmake" -DCMAKE_INSTALL_PREFIX="/home/niklas/projects/uni/coaler" -DCMAKE_POLICY_DEFAULT_CMP0091="NEW" -DCMAKE_BUILD_TYPE="Debug" "/home/niklas/projects/uni/coaler"
```

Copy this part (it will be different for on your PC):

```
-DCMAKE_TOOLCHAIN_FILE="/home/niklas/projects/uni/coaler/build/Debug/generators/conan_toolchain.cmake" -DCMAKE_INSTALL_PREFIX="/home/niklas/projects/uni/coaler" `
```

1. Open your Project in Clion and tell it that it is a CMake project
2. In Clion Goto Settings->Build,Execution...->CMAKE and add the thing you just copied to `CMake Options`.
3. Change "Build directory" from `cmake-build-debug` to `build` (Conan uses this one)
4. You'll get a CMake plugin tab in your plugins -> Run "Reload CMake"
5. If your build worked, this should work without error.

Now you should be able to see the RDKit headers you linked under `External Libraries -> Header Search Paths -> include`
