/*
 *  test_constexpr_verify.cpp - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  test_constexpr_verify is a c++ module
 */

import std;
import lam.linearalgebra;

int main()
{
  // Test constexpr construction and indexing
  constexpr auto test_constexpr = []() {
    lam::linalg::vector<double> v{1.0, 2.0, 3.0};
    return v[1];
  };
  static_assert(test_constexpr() == 2.0);

  // Test constexpr dot product
  constexpr auto test_dot = []() {
    lam::linalg::vector<double> v1{1.0, 2.0};
    lam::linalg::vector<double> v2{3.0, 4.0};
    return v1.dot(v2);
  };
  static_assert(test_dot() == 11.0);

  // Test constexpr norm2
  constexpr auto test_norm2 = []() {
    lam::linalg::vector<double> v{3.0, 4.0};
    return v.norm2();
  };
  static_assert(test_norm2() == 25.0);

  // Test constexpr norm (sqrt)
  constexpr auto test_norm = []() {
    lam::linalg::vector<double> v{3.0, 4.0};
    return v.norm();
  };
  static_assert(test_norm() == 5.0);

  // Matrix accessors, views, transpose, equality and label, in both layouts
  constexpr auto test_matrix_members = []<lam::linalg::storage_layout L>() {
    lam::linalg::matrix<double, L> m(std::vector<std::vector<double>>{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}});
    if (m.row(1)[2] != 6.0 || m.col(2)[0] != 3.0)
      return 1;
    double row_sum = 0.0;
    double col_sum = 0.0;
    for (auto r : m.rows_range())
      for (double x : r)
        row_sum += x;
    for (auto c : m.cols_range())
      for (double x : c)
        col_sum += x;
    if (row_sum != 21.0 || col_sum != 21.0)
      return 2;
    auto t = m.transpose();
    if (t.rows() != 3 || t.cols() != 2 || t[2, 1] != 6.0 || t[0, 1] != 4.0)
      return 3;
    auto copy = m;
    if (!(copy == m) || (copy[0, 0] = 9.0, copy == m))
      return 4;
    m.set_label("m");
    if (m.get_label() != "m")
      return 5;
    return 0;
  };
  static_assert(test_matrix_members.operator()<lam::linalg::storage_layout::row_major>() == 0);
  static_assert(test_matrix_members.operator()<lam::linalg::storage_layout::col_major>() == 0);

  // LU is pure arithmetic -- no BLAS, no LAPACK -- so it belongs at compile time too.
  // Returns the number of the first failing check.
  auto test_lu = [] consteval {
    using mat = lam::linalg::matrix<double, lam::linalg::storage_layout::row_major>;

    // det = 1*(4*6-5*0) - 2*(3*6-5*1) + 0 = 24 - 26 = -2
    mat a(3, 3);
    a[0, 0] = 1.0;
    a[0, 1] = 2.0;
    a[0, 2] = 0.0;
    a[1, 0] = 3.0;
    a[1, 1] = 4.0;
    a[1, 2] = 5.0;
    a[2, 0] = 1.0;
    a[2, 1] = 0.0;
    a[2, 2] = 6.0;

    const double det = lam::linalg::crout_lu_det(a);
    if (det > -2.0 + 1e-12 || det < -2.0 - 1e-12)
      return 1;
    const double ddet = lam::linalg::doolittle_lu_det(a);
    if (ddet > -2.0 + 1e-12 || ddet < -2.0 - 1e-12)
      return 2;

    auto decomp = lam::linalg::crout_lu(a);
    if (!decomp.has_value())
      return 3;

    // A x = b with b = A * (1, 1, 1) must give back (1, 1, 1).
    lam::linalg::vector<double> b({3.0, 12.0, 7.0});
    auto x = lam::linalg::crout_lu_solve(a, b);
    for (std::size_t i = 0; i < 3; ++i)
      if (x[i] > 1.0 + 1e-10 || x[i] < 1.0 - 1e-10)
        return 4;

    // A * A^-1 == I
    auto inv = lam::linalg::crout_lu_inv(a);
    auto prod = a * inv;
    for (std::size_t i = 0; i < 3; ++i)
      for (std::size_t j = 0; j < 3; ++j)
      {
        const double want = (i == j) ? 1.0 : 0.0;
        if (prod[i, j] > want + 1e-10 || prod[i, j] < want - 1e-10)
          return 5;
      }

    // crout_lu substitutes a tiny pivot, so has_value() does NOT mean invertible; only an all-zero row fails.
    mat s(2, 2);
    s[0, 0] = 1.0;
    s[0, 1] = 2.0;
    s[1, 0] = 2.0;
    s[1, 1] = 4.0;
    if (!lam::linalg::crout_lu(s).has_value())
      return 6;
    const double sdet = lam::linalg::crout_lu_det(s);
    if (sdet > 1e-15 || sdet < -1e-15)
      return 7;

    // An all-zero row is the case that does fail.
    mat z(2, 2);
    z[0, 0] = 1.0;
    z[0, 1] = 2.0;
    if (lam::linalg::crout_lu(z).has_value())
      return 8;

    // Doolittle does not pivot, so a zero pivot must return nullopt rather than divide by it.
    // This matrix is perfectly invertible; it just needs a row swap Doolittle will not do.
    mat np(2, 2);
    np[0, 0] = 0.0;
    np[0, 1] = 1.0;
    np[1, 0] = 1.0;
    np[1, 1] = 0.0;
    if (lam::linalg::doolittle_lu(np).has_value())
      return 9;
    if (lam::linalg::doolittle_lu(s).has_value()) // singular: zero pivot at k=1
      return 10;
    // crout_lu pivots, so it handles the same matrix fine: det = -1.
    const double npdet = lam::linalg::crout_lu_det(np);
    if (npdet > -1.0 + 1e-12 || npdet < -1.0 - 1e-12)
      return 11;

    // Scaling cannot make an invertible matrix singular, so the pivot tests must be relative.
    auto scaled = [](double c) consteval {
      mat m(2, 2);
      m[0, 0] = 1.0 * c;
      m[0, 1] = 2.0 * c;
      m[1, 0] = 3.0 * c;
      m[1, 1] = 4.0 * c;
      return m;
    };

    // Entries near epsilon, det = -2e-34: representable, yet an absolute threshold calls it singular.
    if (!lam::linalg::crout_lu(scaled(1e-17)).has_value())
      return 12;
    if (!lam::linalg::doolittle_lu(scaled(1e-17)).has_value())
      return 13;

    // det(c*A) = c^2 * det(A) for a 2x2, checked relatively across twenty orders of magnitude.
    for (double c : {1e-10, 1e-3, 1.0, 1e3, 1e10})
    {
      const double want = -2.0 * c * c;
      const double got = lam::linalg::crout_lu_det(scaled(c));
      const double slack = (want < 0.0 ? -want : want) * 1e-12;
      if (got - want > slack || want - got > slack)
        return 14;
    }

    // The threshold must follow the scalar type; casting through double misfits float at both ends.
    using fmat = lam::linalg::matrix<float, lam::linalg::storage_layout::row_major>;
    fmat f(2, 2);
    f[0, 0] = 1e-20f;
    f[0, 1] = 2e-20f;
    f[1, 0] = 3e-20f;
    f[1, 1] = 4e-20f;
    if (!lam::linalg::crout_lu(f).has_value())
      return 15;

    return 0;
  };
  static_assert(test_lu() == 0);

  std::println("Constexpr tests passed!");
  return 0;
}
