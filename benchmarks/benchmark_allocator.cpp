/*
 *  benchmark_allocator.cpp - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  benchmark_allocator is a c++ module
 */

import std;
import lam.linearalgebra;

using namespace std::chrono;
using namespace lam::linalg;

// Prevent optimization
volatile double sink;

void bench_heap(int iterations)
{
  using vec_t = vector<double>; // Default allocator (heap)
  vec_t v1({1.0, 2.0, 3.0});
  vec_t v2({4.0, 5.0, 6.0});

  auto start = steady_clock::now();
  for (int i = 0; i < iterations; ++i)
  {
    vec_t v3 = v1 + v2;
    sink = v3[0];
  }
  auto end = steady_clock::now();

  auto dur = duration_cast<microseconds>(end - start).count();
  std::println("Heap Allocator: {} us for {} iterations", dur, iterations);
}

void bench_pmr_monotonic(int iterations)
{
  using vec_t = vector<double, std::pmr::polymorphic_allocator<double>>;

  // Large buffer to avoid upstream allocations during benchmark
  // 3 doubles * 8 bytes * 1000000 iterations = ~24MB
  // Using 64MB buffer
  std::vector<std::byte> buffer(64 * 1024 * 1024);
  std::pmr::monotonic_buffer_resource pool{buffer.data(), buffer.size()};
  std::pmr::polymorphic_allocator<double> alloc{&pool};

  vec_t v1({1.0, 2.0, 3.0}, alloc);
  vec_t v2({4.0, 5.0, 6.0}, alloc);

  auto start = steady_clock::now();
  for (int i = 0; i < iterations; ++i)
  {
    // Monotonic buffer just bumps pointer.
    // NOTE: In real loop, you would typically reset() the buffer periodically
    // if it's scratch memory. Here we just consume linear memory to measure allocation speed.
    vec_t v3 = v1 + v2;
    sink = v3[0];
  }
  auto end = steady_clock::now();

  auto dur = duration_cast<microseconds>(end - start).count();
  std::println("PMR Monotonic:  {} us for {} iterations", dur, iterations);
}

void bench_pmr_pool(int iterations)
{
  using vec_t = vector<double, std::pmr::polymorphic_allocator<double>>;
  // unsynchronized_pool_resource is good for frequent alloc/dealloc of same size
  std::pmr::unsynchronized_pool_resource pool;
  std::pmr::polymorphic_allocator<double> alloc{&pool};

  vec_t v1({1.0, 2.0, 3.0}, alloc);
  vec_t v2({4.0, 5.0, 6.0}, alloc);

  auto start = steady_clock::now();
  for (int i = 0; i < iterations; ++i)
  {
    vec_t v3 = v1 + v2;
    sink = v3[0];
  }
  auto end = steady_clock::now();

  auto dur = duration_cast<microseconds>(end - start).count();
  std::println("PMR Pool:       {} us for {} iterations", dur, iterations);
}

int main()
{
  std::println("Benchmarking Vector Addition...");
  int iterations = 1000000;

  bench_heap(iterations);
  bench_pmr_monotonic(iterations);
  bench_pmr_pool(iterations);

  return 0;
}
