/*
 *  benchmark_matrix_mult.cpp - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  benchmark_matrix_mult is a c++ module
 */

import std;
import lam.linearalgebra;

using namespace lam::linalg;
using namespace std::chrono;

template<typename Func>
void measure(const char* name, int N, Func&& f)
{
  // Warmup
  // f();

  auto start = steady_clock::now();
  f();
  auto end = steady_clock::now();
  auto duration = duration_cast<milliseconds>(end - start);

  double seconds = duration.count() / 1000.0;
  // 2 * N^3 FLOPs
  double gflops = (2.0 * std::pow(N, 3)) / (seconds * 1e9);

  std::cout << name << ": " << duration.count() << " ms (" << gflops << " GFLOPS)" << std::endl;
}

int main()
{
  constexpr int N = 1024;
  std::cout << "Benchmarking Matrix Multiplication (N=" << N << ")..." << std::endl;

  // 1. Matrix A (Row Major)
  matrix<double> A(N, N);

  // 2. Matrix B (Row Major) - Baseline
  matrix<double> B(N, N);

  // 3. Matrix C (Col Major) - Target for Optimization
  // Physically, C's data is arranged such that Row*Col should be very fast with the right loop.
  matrix<double, std::allocator<double>, storage_layout::col_major> C(N, N);

  // Fill with dummy data
  for (int i = 0; i < N; ++i)
  {
    for (int j = 0; j < N; ++j)
    {
      A[i, j] = 1.0;
      B[i, j] = 1.0;
      C[i, j] = 1.0;
    }
  }

  // Case 1: Row * Row
  // Uses standard i, k, j loop (inner loop j is seq for B).
  measure("RowMajor * RowMajor (Standard i,k,j)", N, [&]() {
    auto res = A * B;
    // prevent optimization
    if (res.rows() == 0)
      std::abort();
  });

  // Case 2: Row * Col
  // Uses optimized i, j, k loop (inner loop k is seq for A and C).
  // This should be comparable or faster, definitely not slower.
  // If it were using i, k, j, it would be SLOW (strided access on C).
  measure("RowMajor * ColMajor (Context Optimized i,j,k)", N, [&]() {
    auto res = A * C;
    if (res.rows() == 0)
      std::abort();
  });

  // 4. Matrix D (Col Major)
  matrix<double, std::allocator<double>, storage_layout::col_major> D(N, N);
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j)
      D[i, j] = 1.0;

  // Case 3: Col * Col
  // Uses default i, k, j loop.
  // res[i, j] writes RowMajor (seq).
  // D[k, j] reads ColMajor (j varies -> strided). SLOW?
  measure("ColMajor * ColMajor (Context Optimized k,j,i)", N, [&]() {
    auto res = C * D;
    if (res.rows() == 0)
      std::abort();
  });

  // Case 4: Col * Row
  // Uses default i, k, j loop.
  // res[i, j] writes RowMajor (seq).
  // B[k, j] reads RowMajor (seq).
  // A[i, k] reads ColMajor (strided? no k varies -> contiguous).
  // Should be FAST.
  measure("ColMajor * RowMajor (Default i,k,j)", N, [&]() {
    auto res = C * B;
    if (res.rows() == 0)
      std::abort();
  });

  return 0;
}
