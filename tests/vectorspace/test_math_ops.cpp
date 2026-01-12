/*
 *  test_math_ops.cpp - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  test_math_ops is a c++ module
 */

import std;
import lam.linearalgebra;

// Helper for approximate equality
constexpr bool approx_eq(double a, double b, double tol = 1e-9)
{
  double diff = a - b;
  return (diff < 0 ? -diff : diff) < tol;
}

// Simple check function that throws on failure
void check(bool condition, std::string_view msg)
{
  if (!condition)
  {
    throw std::runtime_error(std::string(msg));
  }
}

int main()
{
  using lam::linalg::vector;

  // ===== Scalar Division =====
  {
    vector<double> v{6.0, 9.0, 12.0};
    v /= 3.0;
    check(v[0] == 2.0 && v[1] == 3.0 && v[2] == 4.0, "operator/= failed");
    std::println("operator/= passed");
  }
  {
    vector<double> v{6.0, 9.0, 12.0};
    auto result = v / 3.0;
    check(result[0] == 2.0 && result[1] == 3.0 && result[2] == 4.0, "operator/ failed");
    std::println("operator/ passed");
  }

  // ===== Normalize (in-place) =====
  {
    vector<double> v{3.0, 4.0};
    v.normalize();
    check(approx_eq(v[0], 0.6) && approx_eq(v[1], 0.8), "normalize values failed");
    check(approx_eq(v.norm(), 1.0), "normalize norm failed");
    std::println("normalize() passed");
  }

  // ===== Angle =====
  {
    vector<double> v1{1.0, 0.0};
    vector<double> v2{0.0, 1.0};
    double ang = v1.angle(v2);
    constexpr double pi_over_2 = 1.5707963267948966;
    check(approx_eq(ang, pi_over_2, 1e-6), "angle() failed");
    std::println("angle() passed");
  }
  {
    // Parallel vectors should have angle 0
    vector<double> v1{1.0, 2.0};
    vector<double> v2{2.0, 4.0};
    check(approx_eq(v1.angle(v2), 0.0, 1e-6), "angle() parallel failed");
    std::println("angle() parallel passed");
  }

  // ===== Project =====
  {
    vector<double> v{3.0, 4.0};
    vector<double> onto{1.0, 0.0};
    auto proj = v.project(onto);
    check(approx_eq(proj[0], 3.0) && approx_eq(proj[1], 0.0), "project() failed");
    std::println("project() passed");
  }

  // ===== Reject =====
  {
    vector<double> v{3.0, 4.0};
    vector<double> from{1.0, 0.0};
    auto rej = v.reject(from);
    check(approx_eq(rej[0], 0.0) && approx_eq(rej[1], 4.0), "reject() failed");
    std::println("reject() passed");
  }

  // ===== Reflect =====
  {
    vector<double> v{1.0, -1.0};
    vector<double> normal{0.0, 1.0};
    auto reflected = v.reflect(normal);
    check(approx_eq(reflected[0], 1.0) && approx_eq(reflected[1], 1.0), "reflect() failed");
    std::println("reflect() passed");
  }

  // ===== Lerp =====
  {
    vector<double> a{0.0, 0.0};
    vector<double> b{10.0, 20.0};
    auto mid = a.lerp(b, 0.5);
    check(approx_eq(mid[0], 5.0) && approx_eq(mid[1], 10.0), "lerp() failed");
    std::println("lerp() passed");
  }
  {
    // Free function version
    vector<double> a{0.0, 0.0};
    vector<double> b{10.0, 20.0};
    auto mid = lam::linalg::lerp(a, b, 0.25);
    check(approx_eq(mid[0], 2.5) && approx_eq(mid[1], 5.0), "lerp free function failed");
    std::println("lerp() free function passed");
  }

  // ===== Distance =====
  {
    vector<double> a{0.0, 0.0};
    vector<double> b{3.0, 4.0};
    check(approx_eq(a.distance(b), 5.0), "distance() failed");
    check(approx_eq(lam::linalg::distance(a, b), 5.0), "distance free function failed");
    std::println("distance() passed");
  }

  // ===== is_parallel =====
  {
    vector<double> v1{1.0, 2.0, 3.0};
    vector<double> v2{2.0, 4.0, 6.0};
    vector<double> v3{-1.0, -2.0, -3.0};
    vector<double> v4{1.0, 0.0, 0.0};
    check(v1.is_parallel(v2), "is_parallel positive failed");
    check(v1.is_parallel(v3), "is_parallel anti-parallel failed");
    check(!v1.is_parallel(v4), "is_parallel non-parallel failed");
    std::println("is_parallel() passed");
  }

  // ===== is_orthogonal =====
  {
    vector<double> v1{1.0, 0.0, 0.0};
    vector<double> v2{0.0, 1.0, 0.0};
    vector<double> v3{0.0, 0.0, 1.0};
    vector<double> v4{1.0, 1.0, 0.0};
    check(v1.is_orthogonal(v2), "is_orthogonal xy failed");
    check(v1.is_orthogonal(v3), "is_orthogonal xz failed");
    check(!v1.is_orthogonal(v4), "is_orthogonal non-orthogonal failed");
    std::println("is_orthogonal() passed");
  }

  // ===== Triple Product =====
  {
    vector<double> a{1.0, 0.0, 0.0};
    vector<double> b{0.0, 1.0, 0.0};
    vector<double> c{0.0, 0.0, 1.0};
    // Triple product of unit vectors along axes = 1 (volume of unit cube)
    check(approx_eq(a.triple_product(b, c), 1.0), "triple_product() failed");
    std::println("triple_product() passed");
  }

  // ===== Free Functions =====
  {
    vector<double> a{1.0, 2.0, 3.0};
    vector<double> b{4.0, 5.0, 6.0};
    check(lam::linalg::dot(a, b) == a.dot(b), "dot free function failed");
    std::println("dot() free function passed");
  }
  {
    vector<double> a{1.0, 0.0, 0.0};
    vector<double> b{0.0, 1.0, 0.0};
    auto c = lam::linalg::cross(a, b);
    check(approx_eq(c[0], 0.0) && approx_eq(c[1], 0.0) && approx_eq(c[2], 1.0), "cross free function failed");
    std::println("cross() free function passed");
  }
  {
    vector<double> a{1.0, 0.0};
    vector<double> b{0.0, 1.0};
    constexpr double pi_over_2 = 1.5707963267948966;
    check(approx_eq(lam::linalg::angle(a, b), pi_over_2, 1e-6), "angle free function failed");
    std::println("angle() free function passed");
  }

  // ===== Range Functions =====
  {
    std::vector<double> v1{1.0, 2.0, 3.0};
    std::vector<double> v2{4.0, 5.0, 6.0};
    // dot product: 4 + 10 + 18 = 32
    check(approx_eq(lam::linalg::dot_range(v1, v2), 32.0), "dot_range failed");
    std::println("dot_range passed");
  }
  {
    std::vector<double> v1{0.0, 0.0};
    std::vector<double> v2{3.0, 4.0};
    check(approx_eq(lam::linalg::distance_range(v1, v2), 5.0), "distance_range failed");
    std::println("distance_range passed");
  }


  // ===== Constexpr Tests =====
  constexpr auto test_project = []() {
    vector<double> v{3.0, 4.0};
    vector<double> onto{1.0, 0.0};
    auto proj = v.project(onto);
    return proj[0];
  };
  static_assert(test_project() == 3.0);
  std::println("constexpr project() passed");

  constexpr auto test_lerp = []() {
    vector<double> a{0.0, 0.0};
    vector<double> b{10.0, 20.0};
    auto mid = a.lerp(b, 0.5);
    return mid[0];
  };
  static_assert(test_lerp() == 5.0);
  std::println("constexpr lerp() passed");

  constexpr auto test_distance = []() {
    vector<double> a{0.0, 0.0};
    vector<double> b{3.0, 4.0};
    return a.distance(b);
  };
  static_assert(test_distance() == 5.0);
  std::println("constexpr distance() passed");

  constexpr auto test_range_constexpr = []() {
    std::array<double, 3> a{1.0, 2.0, 3.0};
    std::array<double, 3> b{4.0, 5.0, 6.0};
    return lam::linalg::dot_range(a, b);
  };
  static_assert(test_range_constexpr() == 32.0);
  std::println("constexpr dot_range passed");

  constexpr auto test_distance_range_constexpr = []() {
    std::array<double, 2> a{0.0, 0.0};
    std::array<double, 2> b{3.0, 4.0};
    return lam::linalg::distance_range(a, b);
  };
  static_assert(test_distance_range_constexpr() == 5.0);
  std::println("constexpr distance_range passed");


  std::println("\n✅ All math operations tests passed!");
  return 0;
}
