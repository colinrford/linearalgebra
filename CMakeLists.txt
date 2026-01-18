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

project(lam_linearalgebra VERSION 0.3.260118 LANGUAGES CXX)

# Emscripten needs manual __CMAKE::CXX23 target creation since CMake's
# detection doesn't find the std modules in Emscripten's toolchain
if(EMSCRIPTEN AND NOT TARGET "__CMAKE::CXX23")
  set(_emscripten_libcxx_path "/usr/local/Cellar/emscripten/4.0.23/libexec/llvm/share/libc++/v1")
  if(EXISTS "${_emscripten_libcxx_path}/std.cppm")
    add_library(__cmake_cxx23 STATIC)
    target_sources(__cmake_cxx23 INTERFACE "$<$<STREQUAL:$<TARGET_PROPERTY:TYPE>,STATIC_LIBRARY>:$<TARGET_OBJECTS:__cmake_cxx23>>")
    set_property(TARGET __cmake_cxx23 PROPERTY EXCLUDE_FROM_ALL 1)
    set_property(TARGET __cmake_cxx23 PROPERTY CXX_SCAN_FOR_MODULES 1)
    set_property(TARGET __cmake_cxx23 PROPERTY CXX_MODULE_STD 0)
    target_compile_features(__cmake_cxx23 PUBLIC cxx_std_23)
    target_compile_options(__cmake_cxx23 PRIVATE -Wno-reserved-module-identifier)
    target_include_directories(__cmake_cxx23 PRIVATE
      "${_emscripten_libcxx_path}"
      "/usr/local/Cellar/emscripten/4.0.23/libexec/cache/sysroot/include/c++/v1"
      "/usr/local/Cellar/emscripten/4.0.23/libexec/cache/sysroot/include")
    target_sources(__cmake_cxx23
      PUBLIC
      FILE_SET std TYPE CXX_MODULES
        BASE_DIRS "${_emscripten_libcxx_path}"
        FILES "${_emscripten_libcxx_path}/std.cppm" "${_emscripten_libcxx_path}/std.compat.cppm")
    add_library(__CMAKE::CXX23 ALIAS __cmake_cxx23)
    list(APPEND CMAKE_CXX_COMPILER_IMPORT_STD "23")
    message(STATUS "Manually created __CMAKE::CXX23 target for Emscripten")
  else()
    message(WARNING "Emscripten std.cppm not found at ${_emscripten_libcxx_path}")
  endif()
endif()

# Support centralized cmake modules for in-tree builds
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")
endif()
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(Sanitizers)
include(BlasSetup)

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
  linearalgebra_config.cppm.in
  ${CMAKE_CURRENT_BINARY_DIR}/linearalgebra_config.cppm
)

if(NOT TARGET lam::concepts)
  find_package(lam_concepts REQUIRED)
endif()

set(PROJECT_INCLUDE_DIRS
  ${PROJECT_SOURCE_DIR}/src/vectorspace
  ${PROJECT_SOURCE_DIR}/src/transformations
)

add_library(linearalgebra)
add_library(lam::linearalgebra ALIAS linearalgebra)
set_target_properties(linearalgebra PROPERTIES OUTPUT_NAME "lam_linearalgebra")
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

    src/linalg.cppm
    ${CMAKE_CURRENT_BINARY_DIR}/linearalgebra_config.cppm
)

target_include_directories(linearalgebra PUBLIC "$<BUILD_INTERFACE:${PROJECT_INCLUDE_DIRS}>")
target_link_libraries(linearalgebra PUBLIC lam::concepts)

lam_link_blas(linearalgebra)

add_subdirectory(tests)
add_subdirectory(benchmarks)

install(TARGETS linearalgebra 
  EXPORT linearalgebraTargets
  FILE_SET CXX_MODULES DESTINATION lib/cxx/modules/linearalgebra
)

install(EXPORT linearalgebraTargets
  FILE linearalgebraTargets.cmake
  NAMESPACE lam::
  DESTINATION lib/cmake/lam_linearalgebra
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
  DESTINATION lib/cmake/lam_linearalgebra
)
