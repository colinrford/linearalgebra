# Conan 2.x recipe for lam.linearalgebra.
#
# Local development:
#   conan create . --profile <your-profile>
#
# Consumer projects can then declare:
#   requires = "lam_linearalgebra/<version>"   # version comes from the VERSION file
#
# NOTE: linearalgebra has no LICENSE yet, so this recipe omits the `license`
# field and license packaging. Add both when the project is licensed. CPS
# metadata (install(PACKAGE_INFO)) is likewise deferred — the recipe ships only
# the legacy *Config.cmake the hand-rolled install emits. The recipe delegates
# entirely to CMake.

from conan import ConanFile
from conan.tools.cmake import CMake, CMakeToolchain, cmake_layout
from conan.tools.files import load
import os

class LamLinearalgebraConan(ConanFile):
    name = "lam_linearalgebra"
    author = "Colin Ford"
    url = "https://github.com/colinrford/linearalgebra"
    description = (
        "C++23-module linear algebra (vector spaces, matrices, decompositions, "
        "SVD) with optional BLAS/Accelerate backends for the lam project."
    )
    topics = ("c++23", "modules", "linear-algebra", "matrix", "blas", "lam")

    settings = "os", "compiler", "build_type", "arch"
    package_type = "static-library"

    def set_version(self):
        # Single source of truth: the top-level VERSION file, which CMakeLists
        # also reads. Editing VERSION updates both the CMake project version and
        # this recipe — they never drift.
        self.version = load(
            self, os.path.join(self.recipe_folder, "VERSION")
        ).strip()

    def requirements(self):
        # lam.linearalgebra imports lam.concepts (algebraic/numeric concepts).
        # Mirrors target_link_libraries(... lam_concepts::concepts). Pinned to
        # the most recent concepts tag rather than a float range.
        self.requires("lam_concepts/0.1.260603")

    # BLAS/Accelerate is detected from the system at configure time (BlasSetup
    # .cmake); it is not a Conan dependency. On macOS the Accelerate framework
    # is part of the platform.

    exports_sources = (
        "VERSION",
        "CMakeLists.txt",
        "linearalgebra_config.cppm.in",
        "src/*",
        "cmake/*",
        "README.md",
    )

    def layout(self):
        cmake_layout(self)

    def validate(self):
        cppstd = self.settings.compiler.cppstd
        if cppstd is not None:
            std = int(str(cppstd).replace("gnu", ""))
            if std < 23:
                raise Exception(
                    "lam_linearalgebra requires C++23 (compiler.cppstd >= 23)."
                )

    def generate(self):
        tc = CMakeToolchain(self)
        # Mirror what the project's CMakeLists already assumes.
        tc.cache_variables["CMAKE_CXX_STANDARD"] = "23"
        tc.cache_variables["CMAKE_CXX_SCAN_FOR_MODULES"] = "ON"
        # tests/ benchmarks/ aren't in exports_sources.
        tc.cache_variables["LAM_LINEARALGEBRA_BUILD_TESTS"] = "OFF"
        tc.cache_variables["LAM_LINEARALGEBRA_BUILD_BENCHMARKS"] = "OFF"

        # C++23 module deps: do NOT use a Conan deps generator here. Neither
        # CMakeDeps nor CMakeConfigDeps reproduces an imported target's
        # `FILE_SET CXX_MODULES`, so `import lam.concepts;` can't find the
        # module interface. Instead, point find_package at the dependency's OWN
        # installed config inside its Conan package folder — the real
        # lam_conceptsConfig.cmake carries the module file set. CMakeLists
        # already calls find_package(lam_concepts) via lam_find_dependency.
        prefix_paths = [
            self.dependencies[dep].package_folder
            for dep in ("lam_concepts",)
        ]
        tc.cache_variables["CMAKE_PREFIX_PATH"] = ";".join(prefix_paths)
        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        cmake = CMake(self)
        cmake.install()

    def package_info(self):
        # Match the CMake export names so find_package(lam_linearalgebra) and the
        # Conan-generated CMakeDeps both yield the same target.
        self.cpp_info.set_property("cmake_file_name", "lam_linearalgebra")
        self.cpp_info.set_property(
            "cmake_target_name", "lam_linearalgebra::linearalgebra"
        )
        self.cpp_info.libs = ["lam_linearalgebra"]
        self.cpp_info.builddirs = [
            os.path.join("lib", "cmake", "lam_linearalgebra"),
        ]
