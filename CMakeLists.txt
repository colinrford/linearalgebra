cmake_minimum_required(VERSION 3.31.6)

# C++23 import-std UUID ladder + standard settings (must precede project()).
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/LamCxxModuleStd.cmake)

file(READ "${CMAKE_CURRENT_SOURCE_DIR}/VERSION" _lam_linearalgebra_version)
string(STRIP "${_lam_linearalgebra_version}" _lam_linearalgebra_version)
project(lam_linearalgebra VERSION ${_lam_linearalgebra_version} LANGUAGES CXX)
set_property(DIRECTORY APPEND PROPERTY
  CMAKE_CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/VERSION")

# Emscripten: hand-build the __CMAKE::CXX23 std-modules target (no-op elsewhere).
# Only the vcpp-consumed lam modules (concepts, symbols, linearalgebra) carry
# this helper; must run after project() since it creates a target.
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/LamEmscriptenStd.cmake)
lam_setup_emscripten_std()

# Support centralized cmake modules for in-tree builds
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")
endif()
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(Sanitizers)
include(BlasSetup)
include(LamConstexprLimit)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/LamModule.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/LinalgHelpers.cmake)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Map the BLAS/Accelerate booleans into the "true"/"false" strings the config
# partition template expects. These are passed to lam_add_config_partition via
# EXTRA_DEFINES below (LAM_USE_BLAS_BOOL / LAM_BLAS_BACKEND come from BlasSetup).
macro(map_bool_to_string VAR RESULT)
  if(${VAR})
    set(${RESULT} "true")
  else()
    set(${RESULT} "false")
  endif()
endmacro()

map_bool_to_string(LAM_USE_ACCELERATE LAM_USE_ACCELERATE_BOOL)

# lam.concepts dependency (find_package + the BMI-synthesis property dance).
lam_find_dependency(lam_concepts)

# Camp-A test gating: options default ON for a normal dev build; the conanfile
# flips them OFF since tests/ and benchmarks/ aren't in exports_sources.
option(LAM_LINEARALGEBRA_BUILD_TESTS      "Build the unit tests." ON)
option(LAM_LINEARALGEBRA_BUILD_BENCHMARKS "Build the benchmarks." ON)

set(PROJECT_INCLUDE_DIRS
  ${PROJECT_SOURCE_DIR}/src/vectorspace
  ${PROJECT_SOURCE_DIR}/src/transformations
)

add_library(linearalgebra STATIC)
add_library(lam::linearalgebra ALIAS linearalgebra)
add_library(lam_linearalgebra::linearalgebra ALIAS linearalgebra)
set_target_properties(linearalgebra PROPERTIES
  OUTPUT_NAME "lam_linearalgebra"
  EXPORT_NAME linearalgebra)
target_sources(linearalgebra
  PUBLIC
  FILE_SET CXX_MODULES
  BASE_DIRS ${PROJECT_SOURCE_DIR}/src ${CMAKE_CURRENT_BINARY_DIR}
  FILES
    src/vectorspace/vectorspace.cppm
    src/vectorspace/vectorspace_exceptions.cppm
    src/vectorspace/vectorspace_algorithms.cppm
    src/vectorspace/vectorspace_vector.cppm
    src/vectorspace/vectorspace_r_vector.cppm
    src/vectorspace/vectorspace_fixed_vector.cppm
    src/transformations/matrix.cppm
    src/transformations/matrix_product.cppm
    src/transformations/matrix_decomposition.cppm
    src/transformations/matrix_eigenvalue.cppm
    src/transformations/matrix_svd.cppm

    src/linalg.cppm
)

# lam.linearalgebra:config partition — generated at build time (git_hash
# refreshes on commit/checkout without a reconfigure) and appended to the
# CXX_MODULES file set. EXTRA_DEFINES carries the BLAS/Accelerate fields beyond
# the standard config set.
lam_add_config_partition(linearalgebra NAME linearalgebra
  EXTRA_DEFINES
    -DLAM_USE_ACCELERATE_BOOL=${LAM_USE_ACCELERATE_BOOL}
    -DLAM_ACCELERATE_BACKEND=${LAM_ACCELERATE_BACKEND}
    -DLAM_USE_BLAS_BOOL=${LAM_USE_BLAS_BOOL}
    -DLAM_BLAS_BACKEND=${LAM_BLAS_BACKEND})

target_include_directories(linearalgebra PUBLIC "$<BUILD_INTERFACE:${PROJECT_INCLUDE_DIRS}>")
target_link_libraries(linearalgebra PUBLIC lam_concepts::concepts)

lam_link_blas(linearalgebra)

enable_testing()

# config_dump smoke test + build-config diagnostic (imports lam.linearalgebra).
# Always built (not gated) — it is the partition's import/link smoke test.
lam_add_config_dump(linearalgebra NAME linearalgebra NAMESPACE lam::linalg::config)

if(LAM_LINEARALGEBRA_BUILD_TESTS)
  add_subdirectory(tests)
endif()
if(LAM_LINEARALGEBRA_BUILD_BENCHMARKS)
  add_subdirectory(benchmarks)
endif()

# --- install ----------------------------------------------------------------
# Hand-rolled install (matches bitvector): package name lam_linearalgebra,
# target namespace lam_linearalgebra::. CPS (install(PACKAGE_INFO)) is
# deferred — only the legacy *Config.cmake is emitted for now.
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

install(TARGETS linearalgebra
  EXPORT lam_linearalgebraTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  FILE_SET CXX_MODULES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/lam/linearalgebra
)

install(EXPORT lam_linearalgebraTargets
  FILE lam_linearalgebraTargets.cmake
  NAMESPACE lam_linearalgebra::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/lam_linearalgebra
)

# Legacy Config.cmake: pulls in the lam.concepts dependency, then the targets.
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/lam_linearalgebraConfig.cmake" [[
include(CMakeFindDependencyMacro)
find_dependency(lam_concepts)
include("${CMAKE_CURRENT_LIST_DIR}/lam_linearalgebraTargets.cmake")
]])

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/lam_linearalgebraConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/lam_linearalgebraConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/lam_linearalgebraConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/lam_linearalgebra
)
