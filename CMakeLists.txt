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
project(linearalgebra VERSION 0.2.260108 LANGUAGES CXX)

# Sanitizer Options
option(LAM_ENABLE_ASAN "Enable AddressSanitizer" OFF)
option(LAM_ENABLE_UBSAN "Enable UndefinedBehaviorSanitizer" OFF)
option(LAM_ENABLE_TSAN "Enable ThreadSanitizer (cannot combine with ASan)" OFF)

if(LAM_ENABLE_TSAN AND LAM_ENABLE_ASAN)
  message(FATAL_ERROR "Cannot enable both TSan and ASan simultaneously")
endif()

if(LAM_ENABLE_ASAN)
  add_compile_options(-fsanitize=address -fno-omit-frame-pointer)
  add_link_options(-fsanitize=address)
endif()

if(LAM_ENABLE_UBSAN)
  add_compile_options(-fsanitize=undefined -fno-omit-frame-pointer)
  add_link_options(-fsanitize=undefined)
endif()

if(LAM_ENABLE_TSAN)
  add_compile_options(-fsanitize=thread -fno-omit-frame-pointer)
  add_link_options(-fsanitize=thread)
endif()

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
    src/transformations/matrix_decomposition.cppm
    src/linalg.cppm
    ${CMAKE_CURRENT_BINARY_DIR}/config.cppm
)

target_include_directories(linearalgebra PUBLIC "$<BUILD_INTERFACE:${PROJECT_INCLUDE_DIRS}>")
target_link_libraries(linearalgebra PUBLIC lam::concepts)

# Example/Test executable for initial verification
add_executable(la_test 
  tests/vectorspace/test_pnorm.cpp
)
target_link_libraries(la_test PRIVATE linearalgebra lam::concepts)

add_executable(la_constexpr_test 
  tests/vectorspace/test_constexpr_verify.cpp
)
target_link_libraries(la_constexpr_test PRIVATE linearalgebra lam::concepts)

add_executable(la_concepts_test 
  tests/vectorspace/test_concepts_verify.cpp
)
target_link_libraries(la_concepts_test PRIVATE linearalgebra lam::concepts)

add_executable(la_math_ops_test 
  tests/vectorspace/test_math_ops.cpp
)
target_link_libraries(la_math_ops_test PRIVATE linearalgebra lam::concepts)

add_executable(la_allocator_test 
  tests/vectorspace/test_allocator_propagation.cpp
)
target_link_libraries(la_allocator_test PRIVATE linearalgebra lam::concepts)
add_executable(la_benchmark_allocator 
  benchmarks/benchmark_allocator.cpp
)
target_link_libraries(la_benchmark_allocator PRIVATE linearalgebra lam::concepts)
add_executable(example_constexpr_ring 
  examples/constexpr_ring.cpp
)
target_link_libraries(example_constexpr_ring PRIVATE linearalgebra lam::concepts)

add_executable(la_benchmark_ranges benchmarks/benchmark_ranges.cpp)
target_link_libraries(la_benchmark_ranges PRIVATE linearalgebra lam::concepts)

# Note: The user specified toolchain should be used during configuration:
# cmake -B build -S . -DCMAKE_TOOLCHAIN_FILE=~/cmake-toolchains/homebrew_llvm.cmake


add_executable(bench_compile_vector benchmarks/compile_time_vector.cpp)
target_link_libraries(bench_compile_vector PRIVATE linearalgebra lam::concepts)

add_executable(bench_compile_r_vector benchmarks/compile_time_r_vector.cpp)
target_link_libraries(bench_compile_r_vector PRIVATE linearalgebra lam::concepts)

# Algorithm Benchmarks
add_executable(bench_algorithms benchmarks/benchmark_algorithms.cpp)
target_link_libraries(bench_algorithms PRIVATE linearalgebra lam::concepts)

add_executable(bench_compile_algo_legacy benchmarks/compile_time_legacy.cpp)
target_link_libraries(bench_compile_algo_legacy PRIVATE linearalgebra lam::concepts)

add_executable(bench_compile_algo_range benchmarks/compile_time_range.cpp)
target_link_libraries(bench_compile_algo_range PRIVATE linearalgebra lam::concepts)

add_executable(bench_compile_iota benchmarks/compile_time_iota.cpp)
target_link_libraries(bench_compile_iota PRIVATE linearalgebra lam::concepts)


# Installation Configuration
include(GNUInstallDirs)

# Library Installation
install(TARGETS linearalgebra
  EXPORT linearalgebraTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  FILE_SET CXX_MODULES DESTINATION ${CMAKE_INSTALL_LIBDIR}/c++/v1
)

install(EXPORT linearalgebraTargets
  FILE linearalgebraTargets.cmake
  NAMESPACE lam::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/lam_linearalgebra
)

include(CMakePackageConfigHelpers)

file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/lam_linearalgebraConfig.cmake" [[
include(CMakeFindDependencyMacro)
find_dependency(lam_concepts)
include("${CMAKE_CURRENT_LIST_DIR}/linearalgebraTargets.cmake")
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
