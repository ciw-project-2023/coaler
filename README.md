# coaler

<img src="logo/coaler-removebg.png" width="150">

Core Aligner

## Building locally

To build the project locally, you need to have a valid GCC installation as well as
python (>3.10 with venv support) installed. You can build the project with the following commands:

```bash
./confgiure
make
sudo make install
````

The configure script will create a virtual environment in the `venv` directory and install the conan package manager.
The make targets will use conan to build the project. The binary will be installed to `/usr/local/bin/coaler`. You might
have to add `/usr/local/bin` to your `PATH` variable.

### Building Containerized
If you want to use Podman to run the application, you can use the following command:

```bash
make container
```

This will build an container image with the name `coaler:latest`. You can run the container with the following command:

```bash
./run.sh
```

The script will mount the current directory to your PWD in the container and run the application passing all arguments to it.

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
