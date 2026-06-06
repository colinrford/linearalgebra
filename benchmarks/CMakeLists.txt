# Benchmarks. linalg_add_executable (cmake/LinalgHelpers.cmake) handles the
# lam::linearalgebra + lam::concepts link and C++23. Benchmarks are not ctests.

linalg_add_executable(linalg_benchmark_allocator    SOURCE benchmark_allocator.cpp)
linalg_add_executable(linalg_benchmark_matrix_mult  SOURCE benchmark_matrix_mult.cpp)

# These exercise the Accelerate/BLAS path directly, so only build them when a
# BLAS backend is active and link it explicitly (BLAS).
if(LAM_USE_ACCELERATE)
  linalg_add_executable(linalg_benchmark_accelerate SOURCE benchmark_accelerate.cpp BLAS)
  linalg_add_executable(linalg_benchmark_svd        SOURCE benchmark_svd.cpp        BLAS)
endif()

linalg_add_executable(linalg_benchmark_ranges       SOURCE benchmark_ranges.cpp)
linalg_add_executable(bench_compile_vector      SOURCE compile_time_vector.cpp)
linalg_add_executable(bench_compile_r_vector    SOURCE compile_time_r_vector.cpp)

# Algorithm benchmarks
linalg_add_executable(bench_algorithms          SOURCE benchmark_algorithms.cpp)
linalg_add_executable(bench_compile_algo_legacy SOURCE compile_time_legacy.cpp)
linalg_add_executable(bench_compile_algo_range  SOURCE compile_time_range.cpp)
linalg_add_executable(bench_compile_iota        SOURCE compile_time_iota.cpp)
