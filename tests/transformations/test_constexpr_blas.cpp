/*
 *  test_constexpr_blas.cpp - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  test_constexpr_blas is a c++ module
 */

import std;
import lam.linearalgebra;
import lam.concepts;

using namespace lam::linalg;

// This function attempts to multiply matrices at compile time.
// If BLAS is enabled (and checked unchecked in constexpr), this should fail
// because cblas_dgemm is not constexpr.
constexpr bool test_constexpr_mult()
{
  matrix<double> A(2, 2);
  A[0, 0] = 1.0;
  A[0, 1] = 0.0;
  A[1, 0] = 0.0;
  A[1, 1] = 1.0;

  matrix<double> B(2, 2);
  B[0, 0] = 2.0;
  B[0, 1] = 0.0;
  B[1, 0] = 0.0;
  B[1, 1] = 2.0;

  auto C = A * B;
  return (C[0, 0] == 2.0);
}

void test_runtime_mult()
{
  std::println("Linear Algebra Config:");
  std::println("  Use BLAS: {}", config::use_blas);
  std::println("  Backend:  {}", config::blas_backend);

  std::println("Running runtime multiplication (double)...");
  {
    matrix<double> A(2, 2);
    A[0, 0] = 1.0;
    A[0, 1] = 2.0;
    A[1, 0] = 3.0;
    A[1, 1] = 4.0;

    matrix<double> B(2, 2);
    B[0, 0] = 1.0;
    B[0, 1] = 0.0;
    B[1, 0] = 0.0;
    B[1, 1] = 1.0; // Identity

    auto C = A * B;

    bool correct = (C[0, 0] == 1.0 && C[0, 1] == 2.0 && C[1, 0] == 3.0 && C[1, 1] == 4.0);

    if (correct)
    {
      std::println("Runtime multiplication (double) Passed.");
    }
    else
    {
      std::println("Runtime multiplication (double) Failed!");
      std::exit(1);
    }
  }

  std::println("Running runtime multiplication (float)...");
  {
    matrix<float> A(2, 2);
    A[0, 0] = 1.0f;
    A[0, 1] = 2.0f;
    A[1, 0] = 3.0f;
    A[1, 1] = 4.0f;

    matrix<float> B(2, 2);
    B[0, 0] = 1.0f;
    B[0, 1] = 0.0f;
    B[1, 0] = 0.0f;
    B[1, 1] = 1.0f; // Identity

    auto C = A * B;

    bool correct = (C[0, 0] == 1.0f && C[0, 1] == 2.0f && C[1, 0] == 3.0f && C[1, 1] == 4.0f);

    if (correct)
    {
      std::println("Runtime multiplication (float) Passed.");
    }
    else
    {
      std::println("Runtime multiplication (float) Failed!");
      std::exit(1);
    }
  }
}

int main()
{
  // 1. Force compile-time evaluation (Verified via static_assert)
  static_assert(test_constexpr_mult(), "Constexpr matrix multiplication failed!");
  std::println("Constexpr BLAS test passed.");

  // 2. Perform runtime evaluation (Should use BLAS)
  test_runtime_mult();

  return 0;
}
