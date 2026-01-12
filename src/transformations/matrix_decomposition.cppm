/*
 *  matrix_decomposition.cppm - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  matrix_decomposition is a c++ module
 */

export module lam.linearalgebra:matrix.decomposition;

import std;
import lam.concepts;
import :vectorspace;
import :transformations.concepts;
import :matrix;

namespace lam::linalg
{

// Helper for floating-point comparison
constexpr bool compare(const double a, const double b)
{
  constexpr double epsilon = 4 * std::numeric_limits<double>::epsilon();
  return std::abs(b - a) < epsilon;
}

/**
 * LU Decomposition result structure.
 * Stores the combined L/U matrix, row permutations, and parity for determinant.
 */
export template<typename T, typename Alloc = std::allocator<T>>
struct LUdcmp
{
  matrix<T, Alloc> lu_decomp;    // Combined L\U matrix (L has unit diagonal)
  vector<T, Alloc> permutations; // Row permutation indices
  int parity;                    // +1 or -1 for determinant sign

  LUdcmp(matrix<T, Alloc>&& lu, vector<T, Alloc>&& perms, int p)
    : lu_decomp{std::move(lu)}, permutations{std::move(perms)}, parity{p}
  {}
};

namespace detail
{
auto make_indexing_set = [](std::size_t n) { return std::views::iota(std::size_t{0}, n); };

auto make_indexing_set_starting_at = [](std::size_t index, std::size_t size) {
  return std::views::iota(index, index + size);
};

auto make_indexing_set_reverse = [](std::size_t n) { return std::views::reverse(make_indexing_set(n)); };

template<typename M>
using matrix_scalar_t = typename M::scalar_type;

template<typename M>
using matrix_alloc_t = typename M::allocator_type;
// Should add SFINAE or default if allocator_type is missing, but for now assuming it exists or fixing compilation if it
// breaks. Actually, to be safe for ForeignMatrix, let's use std::allocator if missing? Simpler: Just rely on T and
// default Alloc for our internal matrix.
} // namespace detail

/**
 * Crout's LU decomposition with partial pivoting.
 * Based on Numerical Recipes, 3rd Ed.
 */
export template<typename M>
  requires concepts::experimental::matrix_c_weak<M, typename M::scalar_type>
auto crout_lu(const M& A) -> std::optional<LUdcmp<typename M::scalar_type>>
{
  using T = typename M::scalar_type;
  using Alloc = std::allocator<T>; // Defaulting to std::allocator for result for simplicity

  if (!A.rows() || A.rows() != A.cols())
    return std::nullopt; // Must be square

  std::size_t imax = 0;
  int parity = 1;
  const std::size_t n = A.rows();
  T big;
  T temp;
  const auto square = detail::make_indexing_set(n);
  vector<T, Alloc> perms(n);
  vector<T, Alloc> implicitScalingPerRow(n);

  // Create working copy of matrix (using our implementation)
  matrix<T, Alloc> lu(n, n);
  for (auto i : square)
    for (auto j : square)
      lu[i, j] = A[i, j]; // Assumes M supports [i,j] which concept guarantees

  // Find implicit scaling for each row
  for (auto i : square)
  {
    big = T{0};
    for (auto j : square)
    {
      if ((temp = std::abs(lu[i, j])) > big)
        big = temp;
    }
    if (compare(static_cast<double>(big), 0.0))
      return std::nullopt; // Matrix is singular

    implicitScalingPerRow[i] = T{1} / big;
  }

  // Crout's algorithm with partial pivoting
  for (auto k : square)
  {
    big = T{0};
    auto someIndices = detail::make_indexing_set_starting_at(k, n - k);
    for (auto i : someIndices)
    {
      temp = implicitScalingPerRow[i] * std::abs(lu[i, k]);
      if (temp > big)
      {
        big = temp;
        imax = i;
      }
    }

    if (k != imax)
    {
      // Swap rows
      for (auto j : square)
      {
        temp = lu[imax, j];
        lu[imax, j] = lu[k, j];
        lu[k, j] = temp;
      }
      parity = -parity;
      implicitScalingPerRow[imax] = implicitScalingPerRow[k];
    }

    perms[k] = static_cast<T>(imax);

    if (compare(static_cast<double>(lu[k, k]), 0.0))
    {
      // Matrix is singular; use tiny pivot
      lu[k, k] = std::numeric_limits<T>::epsilon();
    }

    someIndices = detail::make_indexing_set_starting_at(k + 1, n - (k + 1));
    for (auto i : someIndices)
    {
      temp = lu[i, k] /= lu[k, k];
      for (auto j : someIndices)
        lu[i, j] -= temp * lu[k, j];
    }
  }

  return LUdcmp<T, Alloc>(std::move(lu), std::move(perms), parity);
}

/**
 * Compute determinant using Crout LU decomposition.
 */
export template<typename M>
  requires concepts::experimental::matrix_c_weak<M, typename M::scalar_type>
typename M::scalar_type crout_lu_det(const M& A)
{
  using T = typename M::scalar_type;
  auto decomp = crout_lu(A);
  if (!decomp.has_value())
    return T{0};

  const auto& lu = decomp.value();
  T det = static_cast<T>(lu.parity);
  const auto n = A.rows();

  for (std::size_t i = 0; i < n; ++i)
    det *= lu.lu_decomp[i, i];

  return det;
}

/**
 * Solve linear system Ax = b using Crout LU decomposition.
 */
export template<typename M>
  requires concepts::experimental::matrix_c_weak<M, typename M::scalar_type>
auto crout_lu_solve(const M& A, const vector<typename M::scalar_type>& b) -> vector<typename M::scalar_type>
{
  using T = typename M::scalar_type;
  using Alloc = std::allocator<T>;

  auto decomp_opt = crout_lu(A);
  if (!decomp_opt.has_value())
    return vector<T, Alloc>(b.size()); // Return zero vector on failure

  const auto& decomp = decomp_opt.value();
  const auto& lu = decomp.lu_decomp;
  const auto& perms = decomp.permutations;

  std::size_t ii = 0;
  std::size_t ip;
  std::size_t n = b.size();
  T sum = T{0};
  vector<T, Alloc> x(n);

  // Copy b into x
  for (std::size_t rho = 0; rho < n; ++rho)
    x[rho] = b[rho];

  // Forward substitution with unscrambling
  for (std::size_t i = 0; i < n; ++i)
  {
    ip = static_cast<std::size_t>(perms[i]);
    sum = x[ip];
    x[ip] = x[i];
    if (ii != 0)
    {
      for (std::size_t j = ii - 1; j < i; ++j)
        sum -= lu[i, j] * x[j];
    }
    else if (!compare(static_cast<double>(sum), 0.0))
    {
      ii = i + 1;
    }
    x[i] = sum;
  }

  // Back substitution
  for (std::size_t i = n; i-- > 0;)
  {
    sum = x[i];
    for (std::size_t j = i + 1; j < n; ++j)
      sum -= lu[i, j] * x[j];
    x[i] = sum / lu[i, i];
  }

  return x;
}

/**
 * Compute matrix inverse using Crout LU decomposition.
 */
export template<typename M>
  requires concepts::experimental::matrix_c_weak<M, typename M::scalar_type>
auto crout_lu_inv(const M& A) -> matrix<typename M::scalar_type>
{
  using T = typename M::scalar_type;
  using Alloc = std::allocator<T>;

  const std::size_t n = A.rows();

  // Solve for each column of the identity
  matrix<T, Alloc> inv(n, n);
  for (std::size_t j = 0; j < n; ++j)
  {
    // Create column j of identity
    vector<T, Alloc> e(n);
    for (std::size_t i = 0; i < n; ++i)
      e[i] = (i == j) ? T{1} : T{0};

    // Solve A * x = e_j
    vector<T, Alloc> col = crout_lu_solve(A, e);

    // Copy to inverse matrix
    for (std::size_t i = 0; i < n; ++i)
      inv[i, j] = col[i];
  }

  return inv;
}

/**
 * Doolittle's LU decomposition.
 */
export template<typename M>
  requires concepts::experimental::matrix_c_weak<M, typename M::scalar_type>
auto doolittle_lu(const M& A) -> std::optional<LUdcmp<typename M::scalar_type>>
{
  using T = typename M::scalar_type;
  using Alloc = std::allocator<T>;

  if (!A.rows() || A.rows() != A.cols())
    return std::nullopt;

  const std::size_t n = A.rows();
  const auto square = detail::make_indexing_set(n);

  matrix<T, Alloc> lu(n, n);
  for (auto i : square)
    for (auto j : square)
      lu[i, j] = T{0};

  for (auto k : square)
  {
    const auto subsq = detail::make_indexing_set(k + 1);
    const auto subsqc = detail::make_indexing_set_starting_at(k, n - k);
    const auto subsqcpo = detail::make_indexing_set_starting_at(k + 1, n - (k + 1));

    for (auto j : subsqc)
    {
      T sum = T{0};
      for (auto p : subsq)
        sum += lu[k, p] * lu[p, j];
      lu[k, j] = A[k, j] - sum;
    }

    for (auto i : subsqcpo)
    {
      T sum = T{0};
      for (auto p : subsq)
        sum += lu[i, p] * lu[p, k];
      lu[i, k] = (A[i, k] - sum) / lu[k, k];
    }
  }

  vector<T, Alloc> perms(n);
  for (std::size_t i = 0; i < n; ++i)
    perms[i] = static_cast<T>(i); // Identity permutation

  return LUdcmp<T, Alloc>(std::move(lu), std::move(perms), 1);
}

/**
 * Compute determinant using Doolittle LU decomposition.
 */
export template<typename M>
  requires concepts::experimental::matrix_c_weak<M, typename M::scalar_type>
typename M::scalar_type doolittle_lu_det(const M& A)
{
  using T = typename M::scalar_type;
  auto decomp = doolittle_lu(A);
  if (!decomp.has_value())
    return T{0};

  const auto& lu = decomp.value();
  T det = T{1};
  const auto n = A.rows();

  for (std::size_t i = 0; i < n; ++i)
    det *= lu.lu_decomp[i, i];

  return det;
}

} // namespace lam::linalg
