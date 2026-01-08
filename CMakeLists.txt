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
    src/transformations/matrix.cppm
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



# Note: The user specified toolchain should be used during configuration:
# cmake -B build -S . -DCMAKE_TOOLCHAIN_FILE=~/cmake-toolchains/homebrew_llvm.cmake

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
