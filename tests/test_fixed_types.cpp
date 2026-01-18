/*
 *  test_fixed_types.cpp - Tests for fixed_vector and fixed_matrix
 *    see github.com/colinrford/linearalgebra for more info
 *
 *  This test can be compiled standalone to verify the new types work
 *  without requiring modifications to existing module exports.
 *
 *  Compile with (adjust paths as needed):
 *    clang++ -std=c++23 -fmodules ... test_fixed_types.cpp
 */

import std;
import lam.concepts;
import lam.linearalgebra:vectorspace.fixed_vector;
import lam.linearalgebra:fixed_matrix;

using namespace lam::linalg;

// ============================================================================
// fixed_vector tests
// ============================================================================

constexpr bool test_fixed_vector_construction()
{
  // Default construction
  fixed_vector<double, 3> v1;
  if (v1[0] != 0.0 || v1[1] != 0.0 || v1[2] != 0.0) return false;

  // Variadic construction
  fixed_vector<double, 3> v2{1.0, 2.0, 3.0};
  if (v2[0] != 1.0 || v2[1] != 2.0 || v2[2] != 3.0) return false;

  // Named accessors
  if (v2.x() != 1.0 || v2.y() != 2.0 || v2.z() != 3.0) return false;

  // Size
  if (v2.size() != 3) return false;

  return true;
}

constexpr bool test_fixed_vector_arithmetic()
{
  vec3d a{1.0, 2.0, 3.0};
  vec3d b{4.0, 5.0, 6.0};

  // Addition
  auto c = a + b;
  if (c.x() != 5.0 || c.y() != 7.0 || c.z() != 9.0) return false;

  // Subtraction
  auto d = b - a;
  if (d.x() != 3.0 || d.y() != 3.0 || d.z() != 3.0) return false;

  // Scalar multiplication
  auto e = 2.0 * a;
  if (e.x() != 2.0 || e.y() != 4.0 || e.z() != 6.0) return false;

  auto f = a * 3.0;
  if (f.x() != 3.0 || f.y() != 6.0 || f.z() != 9.0) return false;

  // Scalar division
  auto g = b / 2.0;
  if (g.x() != 2.0 || g.y() != 2.5 || g.z() != 3.0) return false;

  // Negation
  auto h = -a;
  if (h.x() != -1.0 || h.y() != -2.0 || h.z() != -3.0) return false;

  // Compound assignment
  vec3d i{1.0, 1.0, 1.0};
  i += a;
  if (i.x() != 2.0 || i.y() != 3.0 || i.z() != 4.0) return false;

  return true;
}

constexpr bool test_fixed_vector_comparison()
{
  vec3d a{1.0, 2.0, 3.0};
  vec3d b{1.0, 2.0, 3.0};
  vec3d c{1.0, 2.0, 4.0};

  if (!(a == b)) return false;
  if (a == c) return false;

  return true;
}

// ============================================================================
// fixed_matrix tests
// ============================================================================

constexpr bool test_fixed_matrix_construction()
{
  // Default (zero)
  mat4d m1;
  if (m1[0, 0] != 0.0) return false;

  // Identity
  auto m2 = mat4d::identity();
  if (m2[0, 0] != 1.0 || m2[1, 1] != 1.0 || m2[2, 2] != 1.0 || m2[3, 3] != 1.0) return false;
  if (m2[0, 1] != 0.0 || m2[1, 0] != 0.0) return false;

  // Dimensions
  if (m2.rows() != 4 || m2.cols() != 4) return false;
  if (!m2.is_square()) return false;

  return true;
}

constexpr bool test_fixed_matrix_arithmetic()
{
  auto a = mat2d::identity();
  auto b = mat2d::identity();

  // Addition
  auto c = a + b;
  if (c[0, 0] != 2.0 || c[1, 1] != 2.0) return false;

  // Scalar multiplication
  auto d = 3.0 * a;
  if (d[0, 0] != 3.0 || d[1, 1] != 3.0) return false;

  // Matrix multiplication (identity * identity = identity)
  auto e = a * b;
  if (e[0, 0] != 1.0 || e[1, 1] != 1.0 || e[0, 1] != 0.0) return false;

  return true;
}

constexpr bool test_matrix_vector_multiply()
{
  // 2x2 identity * vector = same vector
  auto m = mat2d::identity();
  vec2d v{3.0, 4.0};

  auto result = m * v;
  if (result.x() != 3.0 || result.y() != 4.0) return false;

  // Scaling matrix
  fixed_matrix<double, 2, 2> scale_mat{};
  scale_mat[0, 0] = 2.0;
  scale_mat[1, 1] = 3.0;

  auto scaled = scale_mat * v;
  if (scaled.x() != 6.0 || scaled.y() != 12.0) return false;

  return true;
}

constexpr bool test_transpose()
{
  fixed_matrix<double, 2, 3> m{};
  m[0, 0] = 1.0; m[0, 1] = 2.0; m[0, 2] = 3.0;
  m[1, 0] = 4.0; m[1, 1] = 5.0; m[1, 2] = 6.0;

  auto t = m.transpose();

  // t should be 3x2
  if (t.rows() != 3 || t.cols() != 2) return false;

  if (t[0, 0] != 1.0 || t[0, 1] != 4.0) return false;
  if (t[1, 0] != 2.0 || t[1, 1] != 5.0) return false;
  if (t[2, 0] != 3.0 || t[2, 1] != 6.0) return false;

  return true;
}

constexpr bool test_trace()
{
  auto m = mat3d::identity();
  if (m.trace() != 3.0) return false;

  m[0, 0] = 1.0; m[1, 1] = 2.0; m[2, 2] = 3.0;
  if (m.trace() != 6.0) return false;

  return true;
}

constexpr bool test_graphics_transforms()
{
  // Translation
  auto t = translate(vec3d{1.0, 2.0, 3.0});
  vec4d p{0.0, 0.0, 0.0, 1.0};
  auto tp = t * p;
  if (tp.x() != 1.0 || tp.y() != 2.0 || tp.z() != 3.0 || tp.w() != 1.0) return false;

  // Scale
  auto s = scale(vec3d{2.0, 3.0, 4.0});
  vec4d q{1.0, 1.0, 1.0, 1.0};
  auto sq = s * q;
  if (sq.x() != 2.0 || sq.y() != 3.0 || sq.z() != 4.0) return false;

  return true;
}

// ============================================================================
// Compile-time verification
// ============================================================================

static_assert(test_fixed_vector_construction(), "fixed_vector construction failed");
static_assert(test_fixed_vector_arithmetic(), "fixed_vector arithmetic failed");
static_assert(test_fixed_vector_comparison(), "fixed_vector comparison failed");
static_assert(test_fixed_matrix_construction(), "fixed_matrix construction failed");
static_assert(test_fixed_matrix_arithmetic(), "fixed_matrix arithmetic failed");
static_assert(test_matrix_vector_multiply(), "matrix-vector multiply failed");
static_assert(test_transpose(), "transpose failed");
static_assert(test_trace(), "trace failed");
static_assert(test_graphics_transforms(), "graphics transforms failed");

// Concept satisfaction
static_assert(lam::concepts::experimental::vector_c_weak<vec3d>,
              "vec3d must satisfy vector_c_weak");
static_assert(lam::concepts::experimental::vector_c_weak<vec4f>,
              "vec4f must satisfy vector_c_weak");

int main()
{
  std::println("All fixed_vector tests passed at compile time!");
  std::println("All fixed_matrix tests passed at compile time!");

  // Runtime verification (same tests)
  std::println("\nRuntime tests:");

  std::println("  fixed_vector construction: {}",
               test_fixed_vector_construction() ? "PASS" : "FAIL");
  std::println("  fixed_vector arithmetic:   {}",
               test_fixed_vector_arithmetic() ? "PASS" : "FAIL");
  std::println("  fixed_vector comparison:   {}",
               test_fixed_vector_comparison() ? "PASS" : "FAIL");
  std::println("  fixed_matrix construction: {}",
               test_fixed_matrix_construction() ? "PASS" : "FAIL");
  std::println("  fixed_matrix arithmetic:   {}",
               test_fixed_matrix_arithmetic() ? "PASS" : "FAIL");
  std::println("  matrix-vector multiply:    {}",
               test_matrix_vector_multiply() ? "PASS" : "FAIL");
  std::println("  transpose:                 {}",
               test_transpose() ? "PASS" : "FAIL");
  std::println("  trace:                     {}",
               test_trace() ? "PASS" : "FAIL");
  std::println("  graphics transforms:       {}",
               test_graphics_transforms() ? "PASS" : "FAIL");

  // Demo formatting
  std::println("\nFormatting demo:");
  vec3d v{1.5, 2.5, 3.5};
  std::println("  vec3d: {}", v);

  auto m = mat4d::identity();
  std::println("  mat4d identity:\n{}", m);

  auto t = translate(vec3d{10.0, 20.0, 30.0});
  std::println("  translation matrix:\n{}", t);

  return 0;
}
