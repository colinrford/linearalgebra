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

import :matrix;

namespace lam::linalg
{

// std::abs is not constexpr for floating point before C++26, and LU needs it at compile time.
constexpr auto abs_helper(auto x)
{
  using T = decltype(x);
  return x < T{0} ? -x : x;
}

// Negligible = lost in the rounding of its own scale; n is a dimensionless operation count.
template<typename T>
constexpr bool negligible(const T value, const T scale, const std::size_t n)
{ return abs_helper(value) <= static_cast<T>(n) * std::numeric_limits<T>::epsilon() * scale; }

/**
 * LU Decomposition result structure.
 * Stores the combined L/U matrix, row permutations, and parity for determinant.
 */
export template<typename T, typename Alloc = std::allocator<T>>
struct LUdcmp
{
  matrix<T, storage_layout::row_major, Alloc> lu_decomp; // Combined L\U matrix (L has unit diagonal)
  vector<T, Alloc> permutations;                         // Row permutation indices
  int parity;                                            // +1 or -1 for determinant sign

  constexpr LUdcmp(matrix<T, storage_layout::row_major, Alloc>&& lu, vector<T, Alloc>&& perms, int p)
    : lu_decomp{std::move(lu)}, permutations{std::move(perms)}, parity{p}
  {}
};

namespace detail
{
constexpr auto make_indexing_set = [](std::size_t n) { return std::views::iota(std::size_t{0}, n); };

constexpr auto make_indexing_set_starting_at = [](std::size_t index, std::size_t size) {
  return std::views::iota(index, index + size);
};

constexpr auto make_indexing_set_reverse = [](std::size_t n) { return std::views::reverse(make_indexing_set(n)); };
} // namespace detail

/**
 * Crout's LU decomposition with partial pivoting.
 * Based on Numerical Recipes, 3rd Ed.
 */
export template<typename M>
  requires lam::concepts::experimental::matrix_c_weak<M, typename M::scalar_type>
constexpr auto crout_lu(const M& A) -> std::optional<LUdcmp<typename M::scalar_type, detail::result_alloc_t<M>>>
{
  using T = typename M::scalar_type;
  using Alloc = detail::result_alloc_t<M>;
  const Alloc alloc = detail::result_allocator(A);

  if (!A.rows() || A.rows() != A.cols())
    return std::nullopt; // Must be square

  std::size_t imax = 0;
  int parity = 1;
  const std::size_t n = A.rows();
  T big;
  T temp;
  const auto square = detail::make_indexing_set(n);
  vector<T, Alloc> perms(n, alloc);
  vector<T, Alloc> implicitScalingPerRow(n, alloc);

  // Create working copy of matrix (using our implementation)
  matrix<T, storage_layout::row_major, Alloc> lu(n, n, alloc);
  for (auto i : square)
    for (auto j : square)
      lu[i, j] = A[i, j]; // Assumes M supports [i,j] which concept guarantees

  // Find implicit scaling for each row
  for (auto i : square)
  {
    big = T{0};
    for (auto j : square)
    {
      if ((temp = abs_helper(lu[i, j])) > big)
        big = temp;
    }
    // Exact: a row is unusable only when entirely zero; a uniformly tiny row is still invertible.
    if (big == T{0})
      return std::nullopt; // Matrix is singular

    implicitScalingPerRow[i] = T{1} / big;
  }

  // Crout's algorithm with partial pivoting
  for (auto k : square)
  {
    big = T{0};
    // Must be reset: an all-zero remaining column leaves imax stale and swaps factored rows.
    imax = k;
    auto someIndices = detail::make_indexing_set_starting_at(k, n - k);
    for (auto i : someIndices)
    {
      temp = implicitScalingPerRow[i] * abs_helper(lu[i, k]);
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

    // The reciprocal is row k's own scale, the level the pivot must clear.
    const T row_scale = T{1} / implicitScalingPerRow[k];
    if (negligible(lu[k, k], row_scale, n))
    {
      // Matrix is singular; substitute a pivot that is tiny on this row's scale, not absolutely.
      lu[k, k] = std::numeric_limits<T>::epsilon() * row_scale;
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
  requires lam::concepts::experimental::matrix_c_weak<M, typename M::scalar_type>
constexpr typename M::scalar_type crout_lu_det(const M& A)
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
export template<typename M, typename BAlloc>
  requires lam::concepts::experimental::matrix_c_weak<M, typename M::scalar_type>
constexpr auto crout_lu_solve(const M& A, const vector<typename M::scalar_type, BAlloc>& b)
  -> vector<typename M::scalar_type, detail::result_alloc_t<M>>
{
  using T = typename M::scalar_type;
  using Alloc = detail::result_alloc_t<M>;
  const Alloc alloc = detail::result_allocator(A);

  auto decomp_opt = crout_lu(A);
  if (!decomp_opt.has_value())
    return vector<T, Alloc>(b.size(), alloc); // Return zero vector on failure

  const auto& decomp = decomp_opt.value();
  const auto& lu = decomp.lu_decomp;
  const auto& perms = decomp.permutations;

  std::size_t ii = 0;
  std::size_t ip;
  std::size_t n = b.size();
  T sum = T{0};
  vector<T, Alloc> x(n, alloc);

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
    else if (sum != T{0})
    {
      // Records where b's leading zeros end; an optimization, not a singularity test.
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
  requires lam::concepts::experimental::matrix_c_weak<M, typename M::scalar_type>
constexpr auto crout_lu_inv(const M& A)
  -> matrix<typename M::scalar_type, storage_layout::row_major, detail::result_alloc_t<M>>
{
  using T = typename M::scalar_type;
  using Alloc = detail::result_alloc_t<M>;
  const Alloc alloc = detail::result_allocator(A);

  const std::size_t n = A.rows();
  // Solve for each column of the identity
  matrix<T, storage_layout::row_major, Alloc> inv(n, n, alloc);
  for (std::size_t j = 0; j < n; ++j)
  {
    // Create column j of identity
    vector<T, Alloc> e(n, alloc);
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
  requires lam::concepts::experimental::matrix_c_weak<M, typename M::scalar_type>
constexpr auto doolittle_lu(const M& A) -> std::optional<LUdcmp<typename M::scalar_type, detail::result_alloc_t<M>>>
{
  using T = typename M::scalar_type;
  using Alloc = detail::result_alloc_t<M>;
  const Alloc alloc = detail::result_allocator(A);

  if (!A.rows() || A.rows() != A.cols())
    return std::nullopt;

  const std::size_t n = A.rows();
  const auto square = detail::make_indexing_set(n);

  matrix<T, storage_layout::row_major, Alloc> lu(n, n, alloc);
  for (auto i : square)
    for (auto j : square)
      lu[i, j] = T{0};

  for (auto k : square)
  {
    const auto subsq = detail::make_indexing_set(k + 1);
    const auto subsqc = detail::make_indexing_set_starting_at(k, n - k);
    const auto subsqcpo = detail::make_indexing_set_starting_at(k + 1, n - (k + 1));

    T row_scale = T{0};
    for (auto j : subsqc)
    {
      T sum = T{0};
      for (auto p : subsq)
        sum += lu[k, p] * lu[p, j];
      lu[k, j] = A[k, j] - sum;
      if (abs_helper(lu[k, j]) > row_scale)
        row_scale = abs_helper(lu[k, j]);
    }
    // Doolittle does not pivot: a zero diagonal is fatal, and the divide below would yield inf/nan.
    if (row_scale == T{0} || negligible(lu[k, k], row_scale, n))
      return std::nullopt;

    for (auto i : subsqcpo)
    {
      T sum = T{0};
      for (auto p : subsq)
        sum += lu[i, p] * lu[p, k];
      lu[i, k] = (A[i, k] - sum) / lu[k, k];
    }
  }

  vector<T, Alloc> perms(n, alloc);
  for (std::size_t i = 0; i < n; ++i)
    perms[i] = static_cast<T>(i);

  return LUdcmp<T, Alloc>(std::move(lu), std::move(perms), 1);
}

/**
 * Compute determinant using Doolittle LU decomposition.
 */
export template<typename M>
  requires lam::concepts::experimental::matrix_c_weak<M, typename M::scalar_type>
constexpr typename M::scalar_type doolittle_lu_det(const M& A)
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
