from conan import ConanFile
from conan.tools.cmake import CMakeToolchain, CMake, cmake_layout, CMakeDeps
from conan.tools.env import VirtualRunEnv
import hashlib
import os


class alignerRecipe(ConanFile):
    name = "aligner"
    version = "0.0.1"
    package_type = "application"

    # Optional metadata
    settings = "os", "compiler", "build_type", "arch"

    # Sources are located in the same place as this recipe, copy them to the recipe
    exports_sources = "CMakeLists.txt", "src/*"

    def requirements(self):
        with open(os.path.expanduser("~/.conan2/profiles/default"), 'rb') as profile:
            hash = hashlib.md5(profile.read()).hexdigest()
            self.requires("rdkit/0.0.1@ciw/{}".format(hash))
        self.requires("boost/1.83.0")
        self.requires("catch2/2.13.10")
        self.requires("spdlog/1.12.0")


    def layout(self):
        cmake_layout(self)

    def generate(self):
        ev = VirtualRunEnv(self)
        ev.generate()
        deps = CMakeDeps(self)
        deps.generate()
        tc = CMakeToolchain(self)
        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        cmake = CMake(self)
        cmake.install()

    

    
