/*
 *  test_blas_dispatch.cpp - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  test_blas_dispatch is a c++ module
 */

import std;
import lam.linearalgebra;

void check(bool condition, std::string_view msg)
{
  if (!condition)
  {
    std::println("FAILED: {}", msg);
    throw std::runtime_error(std::string(msg));
  }
}

constexpr bool approx_eq(double a, double b, double tol = 1e-9)
{
  double diff = a - b;
  return (diff < 0 ? -diff : diff) < tol;
}

int main()
{
  using lam::linalg::vector;

  std::println("Testing generic type (int)...");
  {
    vector<int> v1{1, 2, 3};
    vector<int> v2{4, 5, 6};
    // dot: 1*4 + 2*5 + 3*6 = 4 + 10 + 18 = 32
    int d = v1.dot(v2);
    check(d == 32, "int dot product failed");
    std::println("int dot passed: {}", d);

    // norm2: 1+4+9 = 14
    int n2 = v1.norm2();
    check(n2 == 14, "int norm2 failed");
    std::println("int norm2 passed: {}", n2);

    // norm: sqrt(14) ~ 3 (integer sqrt usually not defined well for library unless we check return type)
    // The library norm returns sqrt(norm2(v)). For int, sqrt(int) -> double usually in std::sqrt.
    // But let's check what `norm` returns for int vector.
    // `vector_c_weak` defines scalar_type.
    // `norm` returns `V::scalar_type`.
    // If scalar_type is int, `sqrt_helper` returns int?
    // `sqrt_helper` does `g = (g + x/g)/2`. Integer division.
    // sqrt(14) = 3.74. Integer: 3.
    // 14 -> 7 -> (7 + 14/7)/2 = 4 -> (4 + 14/4)/2 = (4+3)/2 = 3.
    // So it should be 3.
    int n = v1.norm();
    check(n == 3, "int norm failed");
    std::println("int norm passed: {}", n);
  }

  std::println("Testing BLAS-compatible type (double)...");
  {
    vector<double> v1{1.0, 2.0, 3.0};
    vector<double> v2{4.0, 5.0, 6.0};
    double d = v1.dot(v2);
    check(approx_eq(d, 32.0), "double dot product failed");
    std::println("double dot passed: {}", d);

    double n2 = v1.norm2();
    check(approx_eq(n2, 14.0), "double norm2 failed");

    double n = v1.norm();
    check(approx_eq(n, std::sqrt(14.0)), "double norm failed");
  }

  std::println("Testing BLAS-compatible type (float)...");
  {
    vector<float> v1{1.0f, 2.0f};
    vector<float> v2{2.0f, 3.0f};
    // 2 + 6 = 8
    float d = v1.dot(v2);
    check(std::abs(d - 8.0f) < 1e-5f, "float dot product failed");
    std::println("float dot passed: {}", d);
  }


  std::println("Testing non-contiguous (if applicable? r_vector is contiguous) - skipping for now.");

  std::println("All dispatch tests passed!");
  return 0;
}
