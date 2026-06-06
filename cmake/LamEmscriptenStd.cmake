# LamEmscriptenStd.cmake
# Emscripten C++23 std-modules support for the lam submodules that are consumed
# by the vcpp web build (concepts, symbols, linearalgebra). Other lam submodules
# do NOT need this and should not vendor it.
#
# Master copy lives at <workspace>/cmake/LamEmscriptenStd.cmake; the vcpp-
# consumed submodules vendor their own copy under <submodule>/cmake/. Keep the
# copies in sync from the master.
#
# Provides:
#   lam_setup_emscripten_std()  — call AFTER project() to hand-build the
#                                 __CMAKE::CXX23 target on Emscripten builds.

include_guard(GLOBAL)

# Emscripten's Homebrew Cellar version. Override with -DLAM_EMSCRIPTEN_VERSION=...
set(LAM_EMSCRIPTEN_VERSION "5.0.0" CACHE STRING
  "Homebrew emscripten Cellar version used to locate libc++ std-module headers.")

# -----------------------------------------------------------------------------
# lam_setup_emscripten_std()
#
# Emscripten's toolchain doesn't expose the C++23 std modules to CMake's
# compiler detection, so under EMSCRIPTEN we hand-build the __CMAKE::CXX23
# target from a vendored std.cppm / std.compat.cppm. No-op on non-Emscripten
# builds (and a no-op if the target already exists).
#
# MUST be called AFTER project() — it creates a library target. Expects the
# vendored module sources at <workspace>/cmake/emscripten-std-modules, found
# relative to the calling submodule as ../cmake/emscripten-std-modules. The
# Cellar version is taken from LAM_EMSCRIPTEN_VERSION (see above).
# -----------------------------------------------------------------------------
function(lam_setup_emscripten_std)
  if(NOT EMSCRIPTEN OR TARGET __CMAKE::CXX23)
    return()
  endif()

  set(_emscripten_std_modules "${CMAKE_CURRENT_SOURCE_DIR}/../cmake/emscripten-std-modules")
  set(_emscripten_sysroot
    "/usr/local/Cellar/emscripten/${LAM_EMSCRIPTEN_VERSION}/libexec/cache/sysroot/include")
  set(_emscripten_libcxx_headers "${_emscripten_sysroot}/c++/v1")

  if(NOT EXISTS "${_emscripten_std_modules}/std.cppm")
    message(WARNING "Emscripten std.cppm not found at ${_emscripten_std_modules}")
    return()
  endif()

  add_library(__cmake_cxx23 STATIC)
  target_sources(__cmake_cxx23 INTERFACE "$<$<STREQUAL:$<TARGET_PROPERTY:TYPE>,STATIC_LIBRARY>:$<TARGET_OBJECTS:__cmake_cxx23>>")
  set_property(TARGET __cmake_cxx23 PROPERTY EXCLUDE_FROM_ALL 1)
  set_property(TARGET __cmake_cxx23 PROPERTY CXX_SCAN_FOR_MODULES 1)
  set_property(TARGET __cmake_cxx23 PROPERTY CXX_MODULE_STD 0)
  target_compile_features(__cmake_cxx23 PUBLIC cxx_std_23)
  target_compile_options(__cmake_cxx23 PRIVATE -Wno-reserved-module-identifier)
  target_include_directories(__cmake_cxx23 PRIVATE
    "${_emscripten_libcxx_headers}"
    "${_emscripten_sysroot}")
  target_sources(__cmake_cxx23
    PUBLIC
    FILE_SET std TYPE CXX_MODULES
      BASE_DIRS "${_emscripten_std_modules}"
      FILES "${_emscripten_std_modules}/std.cppm" "${_emscripten_std_modules}/std.compat.cppm")
  add_library(__CMAKE::CXX23 ALIAS __cmake_cxx23)

  # Signal import-std availability to the caller's directory scope (the original
  # inline block set this at directory level; PARENT_SCOPE preserves that).
  list(APPEND CMAKE_CXX_COMPILER_IMPORT_STD "23")
  set(CMAKE_CXX_COMPILER_IMPORT_STD "${CMAKE_CXX_COMPILER_IMPORT_STD}" PARENT_SCOPE)
  message(STATUS "Manually created __CMAKE::CXX23 target for Emscripten")
endfunction()
