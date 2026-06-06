# Helper functions for the linearalgebra module's dev targets (tests, benchmarks).
#
# Each of those subdirs builds many small executables that consume the
# lam.linearalgebra C++23 module the same way (link lam::linearalgebra +
# lam::concepts, request C++23). linalg_add_executable() centralizes that
# boilerplate so it lives in one place instead of being re-typed per target.
#
# Modeled on bitvector's bitvector_add_executable / polynomial_nttp's helper.

include_guard(GLOBAL)

# linalg_add_executable(<name>
#   [SOURCE          <file>]    # source file; defaults to <name>.cpp
#   [LIBS            <libs...>] # extra link libraries
#   [DEFINES         <defs...>] # extra PRIVATE compile definitions
#   [BLAS]                      # also call lam_link_blas(<name>)
#   [CONSTEXPR_LIMIT <n>])      # bump the constexpr step/ops budget (heavy CT work)
#
# Always links lam_linearalgebra::linearalgebra + lam_concepts::concepts and
# requests C++23 (the lam_<module>::<module> namespace convention). The global
# CMAKE_CXX_MODULE_STD / CMAKE_CXX_SCAN_FOR_MODULES settings (set in the root
# CMakeLists) supply the module machinery; consuming lam::linearalgebra carries
# its usage requirements.
function(linalg_add_executable NAME)
  cmake_parse_arguments(ARG "BLAS" "SOURCE;CONSTEXPR_LIMIT" "LIBS;DEFINES" ${ARGN})

  if(NOT ARG_SOURCE)
    set(ARG_SOURCE "${NAME}.cpp")
  endif()

  add_executable(${NAME} ${ARG_SOURCE})
  target_link_libraries(${NAME} PRIVATE lam_linearalgebra::linearalgebra lam_concepts::concepts ${ARG_LIBS})
  target_compile_features(${NAME} PRIVATE cxx_std_23)

  if(ARG_DEFINES)
    target_compile_definitions(${NAME} PRIVATE ${ARG_DEFINES})
  endif()
  if(ARG_BLAS)
    lam_link_blas(${NAME})
  endif()
  if(ARG_CONSTEXPR_LIMIT)
    lam_target_constexpr_limit(${NAME} ${ARG_CONSTEXPR_LIMIT})
  endif()
endfunction()

# linalg_add_test(<name> [args as linalg_add_executable])
# Same as linalg_add_executable but also registers the target as a ctest.
function(linalg_add_test NAME)
  linalg_add_executable(${NAME} ${ARGN})
  add_test(NAME ${NAME} COMMAND ${NAME})
endfunction()
