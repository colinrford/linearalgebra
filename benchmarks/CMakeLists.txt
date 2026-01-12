
add_executable(la_benchmark_allocator benchmark_allocator.cpp)
target_link_libraries(la_benchmark_allocator PRIVATE linearalgebra lam::concepts)

add_executable(la_benchmark_matrix_mult benchmark_matrix_mult.cpp)
target_link_libraries(la_benchmark_matrix_mult PRIVATE linearalgebra lam::concepts)

if(APPLE)
  add_executable(la_benchmark_accelerate benchmark_accelerate.cpp)
  target_link_libraries(la_benchmark_accelerate PRIVATE linearalgebra lam::concepts)
  lam_link_blas(la_benchmark_accelerate) # Ensure it links to BLAS too if needed
endif()

# Benchmarks
add_executable(la_benchmark_ranges benchmark_ranges.cpp)
target_link_libraries(la_benchmark_ranges PRIVATE linearalgebra lam::concepts)

add_executable(bench_compile_vector compile_time_vector.cpp)
target_link_libraries(bench_compile_vector PRIVATE linearalgebra lam::concepts)

add_executable(bench_compile_r_vector compile_time_r_vector.cpp)
target_link_libraries(bench_compile_r_vector PRIVATE linearalgebra lam::concepts)

# Algorithm Benchmarks
add_executable(bench_algorithms benchmark_algorithms.cpp)
target_link_libraries(bench_algorithms PRIVATE linearalgebra lam::concepts)

add_executable(bench_compile_algo_legacy compile_time_legacy.cpp)
target_link_libraries(bench_compile_algo_legacy PRIVATE linearalgebra lam::concepts)

add_executable(bench_compile_algo_range compile_time_range.cpp)
target_link_libraries(bench_compile_algo_range PRIVATE linearalgebra lam::concepts)

add_executable(bench_compile_iota compile_time_iota.cpp)
target_link_libraries(bench_compile_iota PRIVATE linearalgebra lam::concepts)
