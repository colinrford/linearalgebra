cmake_minimum_required(VERSION 3.31.6)

if(CMAKE_VERSION VERSION_GREATER_EQUAL 4.0.3)
  set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD "d0edc3af-4c50-42ea-a356-e2862fe7a444")
elseif(CMAKE_VERSION VERSION_GREATER_EQUAL 4.0.0)
  set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD "a9e1cf81-9932-4810-974b-6eccaf14e457")
elseif(CMAKE_VERSION VERSION_GREATER_EQUAL 3.31.6)
  set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD "0e5b6991-d74f-4b3d-a41c-cf096e0b2508")
endif()

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_MODULE_STD ON)
set(CMAKE_CXX_SCAN_FOR_MODULES ON)

project(linearalgebra VERSION 0.3.260110 LANGUAGES CXX)

# Add cmake/ directory to module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Include Modular Configuration
include(Sanitizers)
include(BlasSetup)

# Output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

find_package(Git)
if(GIT_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
else()
  set(GIT_HASH "unknown")
endif()

if(CMAKE_BUILD_TYPE AND CMAKE_BUILD_TYPE MATCHES "^[Dd][Ee][Bb][Uu][Gg]$")
  set(IS_DEBUG "true")
else()
  set(IS_DEBUG "false")
endif()

configure_file(
  config.cppm.in
  ${CMAKE_CURRENT_BINARY_DIR}/config.cppm
)

find_package(lam_concepts REQUIRED)

# Include directories
set(PROJECT_INCLUDE_DIRS
  ${PROJECT_SOURCE_DIR}/src/vectorspace
  ${PROJECT_SOURCE_DIR}/src/transformations
)

# Create library target using C++23 modules
add_library(linearalgebra)
add_library(lam::linearalgebra ALIAS linearalgebra)
target_sources(linearalgebra
  PUBLIC
  FILE_SET CXX_MODULES
  BASE_DIRS ${PROJECT_SOURCE_DIR}/src ${CMAKE_CURRENT_BINARY_DIR}
  FILES
    src/vectorspace/vectorspace.cppm
    src/vectorspace/vectorspace_concepts.cppm
    src/vectorspace/vectorspace_exceptions.cppm
    src/vectorspace/vectorspace_algorithms.cppm
    src/vectorspace/vectorspace_vector.cppm
    src/vectorspace/vectorspace_r_vector.cppm
    src/transformations/transformations_concepts.cppm
    src/transformations/matrix.cppm
    src/transformations/matrix_product.cppm
    src/transformations/matrix_decomposition.cppm

    src/linalg.cppm
    ${CMAKE_CURRENT_BINARY_DIR}/config.cppm
)

target_include_directories(linearalgebra PUBLIC "$<BUILD_INTERFACE:${PROJECT_INCLUDE_DIRS}>")
target_link_libraries(linearalgebra PUBLIC lam::concepts)

# Link appropriate BLAS backend to the library
lam_link_blas(linearalgebra)

# Tests
add_subdirectory(tests)

# Benchmarks
add_subdirectory(benchmarks)

# Installation
install(TARGETS linearalgebra 
  EXPORT linearalgebraTargets
  FILE_SET CXX_MODULES DESTINATION lib/cxx/modules/linearalgebra
)

install(EXPORT linearalgebraTargets
  FILE linearalgebraTargets.cmake
  NAMESPACE lam::
  DESTINATION lib/cmake/linearalgebra
)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "linearalgebraConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

install(FILES "linearalgebraConfigVersion.cmake"
  DESTINATION lib/cmake/linearalgebra
)
