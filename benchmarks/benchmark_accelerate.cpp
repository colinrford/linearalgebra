/*
 *  benchmark_accelerate.cpp - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  benchmark_accelerate is a c++ module
 */

#include <Accelerate/Accelerate.h>
#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>

import lam.linearalgebra;

using namespace lam::linalg;
using namespace std::chrono;

template<typename Func>
void measure(const char* name, int N, Func&& f)
{
  auto start = steady_clock::now();
  f();
  auto end = steady_clock::now();
  auto duration = duration_cast<milliseconds>(end - start);

  double seconds = duration.count() / 1000.0;
  double gflops = (2.0 * std::pow(N, 3)) / (seconds * 1e9);

  std::cout << name << ": " << duration.count() << " ms (" << gflops << " GFLOPS)" << std::endl;
}

int main()
{
  constexpr int N = 1024;
  std::cout << "Benchmarking vs Apple Accelerate (N=" << N << ")..." << std::endl;

  // 1. Setup Our Matrix
  matrix<double> A(N, N);
  matrix<double> B(N, N);
  for (int i = 0; i < N; ++i)
  {
    for (int j = 0; j < N; ++j)
    {
      A[i, j] = 1.0;
      B[i, j] = 1.0;
    }
  }

  // 2. Measure Ours (Baseline)
  measure("lam::matrix (Row*Row)", N, [&]() {
    auto res = A * B;
    if (res.rows() == 0)
      std::abort();
  });

  // 3. Setup Accelerate (CBLAS)
  // cblas_dgemm expects raw pointers.
  // We can use std::vector for raw storage simulation or just allocate arrays.
  std::vector<double> bufA(N * N, 1.0);
  std::vector<double> bufB(N * N, 1.0);
  std::vector<double> bufC(N * N, 0.0);

  // 4. Measure Accelerate
  measure("Apple Accelerate (cblas_dgemm)", N, [&]() {
    // C = alpha*A*B + beta*C
    // Row-major, No Transpose, No Transpose
    cblas_dgemm(::CblasRowMajor, ::CblasNoTrans, ::CblasNoTrans, N, N, N, 1.0, bufA.data(), N, bufB.data(), N, 0.0,
                bufC.data(), N);
  });

  return 0;
}
