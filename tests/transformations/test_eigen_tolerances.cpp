/*
 *  test_eigen_tolerances.cpp - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  test_eigen_tolerances is a c++ module
 */

import std;
import lam.linearalgebra;

namespace linalg = lam::linalg;
using linalg::storage_layout;

using cmat = linalg::matrix<std::complex<double>, storage_layout::row_major>;

int failures = 0;

void check(bool ok, std::string_view what)
{
  if (!ok)
  {
    std::println("FAILED: {}", what);
    ++failures;
  }
}

// Scaling preserves Hermitian-ness, so the test that picks the routine must not be absolute.
void hermitian_detection_is_scale_free()
{
  auto hermitian = [](double c) {
    cmat m(2, 2);
    m[0, 0] = {2.0 * c, 0.0};
    m[0, 1] = {0.0, 1.0 * c};
    m[1, 0] = {0.0, -1.0 * c};
    m[1, 1] = {3.0 * c, 0.0};
    return m;
  };
  // Off-diagonals differ in the 3rd significant digit: unambiguously NOT Hermitian at any scale.
  auto skewed = [](double c) {
    cmat m(2, 2);
    m[0, 0] = {2.0 * c, 0.0};
    m[0, 1] = {0.0, 1.0 * c};
    m[1, 0] = {0.0, -1.002 * c};
    m[1, 1] = {3.0 * c, 0.0};
    return m;
  };

  cmat id(2, 2);
  id[0, 0] = id[1, 1] = {1.0, 0.0};
  id[0, 1] = id[1, 0] = {0.0, 0.0};

  // At 1e-12 the old absolute rule sent the skewed matrix to zhegvd_, which returned a wrong answer.
  for (double c : {1e-12, 1e-6, 1.0, 1e6})
  {
    auto rh = linalg::solve_gep(hermitian(c), id);
    check(rh.success, std::format("hermitian gep failed at scale {:g}", c));
    if (rh.success)
    {
      // Eigenvalues of [[2,i],[-i,3]] are (5 ± sqrt(5))/2, scaled by c.
      const double lo = c * (5.0 - std::sqrt(5.0)) / 2.0;
      const double hi = c * (5.0 + std::sqrt(5.0)) / 2.0;
      check(std::abs(rh.eigenvalues[0].real() - lo) < 1e-9 * c && std::abs(rh.eigenvalues[1].real() - hi) < 1e-9 * c,
            std::format("hermitian gep eigenvalues wrong at scale {:g}", c));
    }

    auto rs = linalg::solve_gep(skewed(c), id);
    check(rs.success, std::format("non-hermitian gep failed at scale {:g}", c));
    if (rs.success)
    {
      // det = 2*3 - (i)(-1.002i) = 4.998, trace = 5, so λ = (5 ± sqrt(25 - 4*4.998))/2.
      const double disc = std::sqrt(25.0 - 4.0 * 4.998);
      const double lo = c * (5.0 - disc) / 2.0;
      const double hi = c * (5.0 + disc) / 2.0;
      check(std::abs(rs.eigenvalues[0].real() - lo) < 1e-9 * c && std::abs(rs.eigenvalues[1].real() - hi) < 1e-9 * c,
            std::format("non-hermitian gep took the Hermitian branch at scale {:g}: got {:g}, {:g} want {:g}, {:g}", c,
                        rs.eigenvalues[0].real(), rs.eigenvalues[1].real(), lo, hi));
    }
  }
}

// Negligible beta is not finite; both negligible is not determined. Neither may depend on scale.
void infinite_eigenvalues_are_classified()
{
  // B singular in its second entry: the pencil has one finite and one infinite eigenvalue.
  auto pencil = [](double c) {
    cmat a(2, 2);
    a[0, 0] = {1.0 * c, 0.0};
    a[0, 1] = a[1, 0] = {0.0, 0.0};
    a[1, 1] = {1.0 * c, 0.0};
    cmat b(2, 2);
    b[0, 0] = {1.0 * c, 0.0};
    b[0, 1] = b[1, 0] = {0.0, 0.0};
    b[1, 1] = {0.0, 0.0};
    return std::pair{a, b};
  };

  // zggev_ does not normalize: at 1e-15 an absolute 1e-14 threshold calls the finite value infinite.
  for (double c : {1e-15, 1e-8, 1.0, 1e8})
  {
    auto [a, b] = pencil(c);
    auto r = linalg::solve_gep(a, b);
    check(r.success, std::format("singular-B gep failed at scale {:g}", c));
    if (!r.success)
      continue;

    check(r.kinds.size() == 2, std::format("kinds not populated at scale {:g}", c));
    if (r.kinds.size() != 2)
      continue;

    std::size_t n_finite = 0;
    std::size_t n_infinite = 0;
    for (std::size_t i = 0; i < 2; ++i)
    {
      if (r.kinds[i] == linalg::eigenvalue_kind::finite)
        ++n_finite;
      else if (r.kinds[i] == linalg::eigenvalue_kind::infinite)
        ++n_infinite;
    }
    check(n_finite == 1 && n_infinite == 1,
          std::format("scale {:g}: expected 1 finite + 1 infinite, got {} + {}", c, n_finite, n_infinite));

    // The finite one is alpha/beta = 1 at every scale.
    for (std::size_t i = 0; i < 2; ++i)
      if (r.kinds[i] == linalg::eigenvalue_kind::finite)
        check(std::abs(r.eigenvalues[i].real() - 1.0) < 1e-9,
              std::format("scale {:g}: finite eigenvalue should be 1, got {:g}", c, r.eigenvalues[i].real()));

    // An entry reported infinite must actually be infinite, and a finite one must not be.
    for (std::size_t i = 0; i < 2; ++i)
    {
      const bool is_inf = std::isinf(r.eigenvalues[i].real());
      check(is_inf == (r.kinds[i] == linalg::eigenvalue_kind::infinite),
            std::format("scale {:g}: kind disagrees with the value at {}", c, i));
    }
  }
}

// Both singular in the same direction gives 0/0: a third outcome, not a default.
void zero_over_zero_is_indeterminate()
{
  for (double c : {1e-8, 1.0, 1e8})
  {
    cmat a(2, 2);
    a[0, 0] = {1.0 * c, 0.0};
    a[0, 1] = a[1, 0] = {0.0, 0.0};
    a[1, 1] = {0.0, 0.0};
    cmat b(2, 2);
    b[0, 0] = {1.0 * c, 0.0};
    b[0, 1] = b[1, 0] = {0.0, 0.0};
    b[1, 1] = {0.0, 0.0};

    auto r = linalg::solve_gep(a, b);
    check(r.success, std::format("0/0 pencil failed at scale {:g}", c));
    if (!r.success || r.kinds.size() != 2)
      continue;

    std::size_t n_indeterminate = 0;
    for (std::size_t i = 0; i < 2; ++i)
      if (r.kinds[i] == linalg::eigenvalue_kind::indeterminate)
        ++n_indeterminate;
    check(n_indeterminate == 1,
          std::format("scale {:g}: expected exactly 1 indeterminate, got {}", c, n_indeterminate));

    // An indeterminate entry must not be reported as a usable number.
    for (std::size_t i = 0; i < 2; ++i)
      if (r.kinds[i] == linalg::eigenvalue_kind::indeterminate)
        check(std::isnan(r.eigenvalues[i].real()),
              std::format("scale {:g}: indeterminate eigenvalue should not be a value", c));
  }
}

// A complex diagonal is not Hermitian and must go to the general routine, whose answer is complex.
void complex_diagonal_is_not_hermitian()
{
  cmat a(2, 2);
  a[0, 0] = {2.0, 0.5};
  a[0, 1] = {0.0, 1.0};
  a[1, 0] = {0.0, -1.0};
  a[1, 1] = {3.0, 0.0};
  cmat id(2, 2);
  id[0, 0] = id[1, 1] = {1.0, 0.0};
  id[0, 1] = id[1, 0] = {0.0, 0.0};

  auto r = linalg::solve_gep(a, id);
  check(r.success, "complex-diagonal gep failed");
  if (!r.success)
    return;
  const bool any_complex = std::abs(r.eigenvalues[0].imag()) > 1e-9 || std::abs(r.eigenvalues[1].imag()) > 1e-9;
  check(any_complex, "a complex diagonal was treated as Hermitian: eigenvalues came back real");
}

// A well-conditioned ordinary problem must be classified entirely finite, at any scale.
void ordinary_pencils_are_all_finite()
{
  for (double c : {1e-8, 1.0, 1e8})
  {
    cmat a(2, 2);
    a[0, 0] = {2.0 * c, 0.0};
    a[0, 1] = {0.0, 1.0 * c};
    a[1, 0] = {0.0, -1.0 * c};
    a[1, 1] = {3.0 * c, 0.0};
    cmat b(2, 2);
    b[0, 0] = b[1, 1] = {1.0, 0.0};
    b[0, 1] = b[1, 0] = {0.0, 0.0};

    auto r = linalg::solve_gep(a, b);
    check(r.success, std::format("ordinary gep failed at scale {:g}", c));
    if (!r.success || r.kinds.size() != 2)
      continue;
    check(r.kinds[0] == linalg::eigenvalue_kind::finite && r.kinds[1] == linalg::eigenvalue_kind::finite,
          std::format("ordinary gep should be all-finite at scale {:g}", c));
  }
}

int main()
{
  hermitian_detection_is_scale_free();
  infinite_eigenvalues_are_classified();
  zero_over_zero_is_indeterminate();
  complex_diagonal_is_not_hermitian();
  ordinary_pencils_are_all_finite();

  if (failures != 0)
  {
    std::println("{} failure(s)", failures);
    return 1;
  }
  std::println("\n✅ eigen tolerance tests passed");
  return 0;
}
