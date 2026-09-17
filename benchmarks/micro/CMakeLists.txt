# benchmarks/micro -- google/benchmark throughput suite (linalg_gb_* targets).
# Reached only when the parent found google/benchmark (find_package QUIET); the
# gb dependency is quarantined here and never touches the library/install surface.

# Benchmark binaries land in bin/benchmarks/, not the flat bin/.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/benchmarks")

# Dependency-clean support MODULE: std + lam.linearalgebra only, NO google/benchmark.
add_library(linalg_bench_support STATIC)
target_sources(linalg_bench_support
  PUBLIC FILE_SET CXX_MODULES FILES bench_support.cppm)
target_link_libraries(linalg_bench_support
  PUBLIC lam_linearalgebra::linearalgebra lam_concepts::concepts)
target_compile_features(linalg_bench_support PUBLIC cxx_std_23)

# Each micro benchmark: existing helper (lam link + C++23) + support module + gb.
# linalg_gb_ prefix flags "needs google/benchmark, slide-grade suite".
linalg_add_executable(linalg_gb_smoke
  SOURCE bm_smoke.cpp
  LIBS linalg_bench_support benchmark::benchmark benchmark::benchmark_main)

linalg_add_executable(linalg_gb_allocator_profiles
  SOURCE bm_allocator_profiles.cpp
  LIBS linalg_bench_support benchmark::benchmark benchmark::benchmark_main)

# Real suites land here next:
# linalg_add_executable(linalg_gb_matrix_mult        SOURCE bm_matrix_mult.cpp ...)
# linalg_add_executable(linalg_gb_vector_ops         SOURCE bm_vector_ops.cpp ...)
