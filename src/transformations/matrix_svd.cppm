/*
 *  matrix_svd.cppm - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  matrix_svd is a c++ module
 */

export module lam.linearalgebra:matrix.svd;

import std;
import lam.concepts;
import :config;
import :vectorspace;
import :matrix;
import :matrix.product;

namespace lam::linalg
{

#ifdef LAM_USE_BLAS
extern "C"
{
  void sgesvd_(const char* jobu, const char* jobvt, const int* m, const int* n, float* a, const int* lda, float* s,
               float* u, const int* ldu, float* vt, const int* ldvt, float* work, const int* lwork, int* info);

  void dgesvd_(const char* jobu, const char* jobvt, const int* m, const int* n, double* a, const int* lda, double* s,
               double* u, const int* ldu, double* vt, const int* ldvt, double* work, const int* lwork, int* info);

  void cgesvd_(const char* jobu, const char* jobvt, const int* m, const int* n, void* a, const int* lda, float* s,
               void* u, const int* ldu, void* vt, const int* ldvt, void* work, const int* lwork, float* rwork,
               int* info);

  void zgesvd_(const char* jobu, const char* jobvt, const int* m, const int* n, void* a, const int* lda, double* s,
               void* u, const int* ldu, void* vt, const int* ldvt, void* work, const int* lwork, double* rwork,
               int* info);
}
#endif

// NOLINTBEGIN(misc-non-private-member-variables-in-classes) — intentional public result struct
export template<typename T, typename Alloc = std::allocator<T>>
struct svd_result
{
  matrix<T, Alloc> u;
  vector<std::conditional_t<std::is_same_v<T, std::complex<double>>, double,
                            std::conditional_t<std::is_same_v<T, std::complex<float>>, float, T>>,
         typename std::allocator_traits<Alloc>::template rebind_alloc<
           std::conditional_t<std::is_same_v<T, std::complex<double>>, double,
                              std::conditional_t<std::is_same_v<T, std::complex<float>>, float, T>>>>
    s;
  matrix<T, Alloc> vt;

  svd_result(matrix<T, Alloc>&& u_in, decltype(s)&& s_in, matrix<T, Alloc>&& vt_in)
    : u(std::move(u_in)), s(std::move(s_in)), vt(std::move(vt_in))
  {}
}; // NOLINTEND(misc-non-private-member-variables-in-classes)

export template<typename M>
  requires lam::concepts::experimental::matrix_c_weak<M, typename M::scalar_type>
auto svd(const M& A) -> std::optional<svd_result<typename M::scalar_type>>
{
  using T = typename M::scalar_type;
  using Alloc = std::allocator<T>;
  using RealT = std::conditional_t<std::is_same_v<T, std::complex<double>>, double,
                                   std::conditional_t<std::is_same_v<T, std::complex<float>>, float, T>>;
  using RealAlloc = typename std::allocator_traits<Alloc>::template rebind_alloc<RealT>;

#ifdef LAM_USE_BLAS
  if constexpr (std::is_same_v<T, float> || std::is_same_v<T, double> || std::is_same_v<T, std::complex<float>> ||
                std::is_same_v<T, std::complex<double>>)
  {
    int m = static_cast<int>(A.rows());
    int n = static_cast<int>(A.cols());
    if (m == 0 || n == 0)
    {
      return std::nullopt;
    }

    matrix<T, Alloc> a_copy(m, n);
    for (int i = 0; i < m; ++i)
    {
      for (int j = 0; j < n; ++j)
      {
        a_copy[i, j] = A[i, j];
      }
    }

    int lap_m = n;
    int lap_n = m;
    int lda = n;

    matrix<T, Alloc> u(m, m);
    matrix<T, Alloc> vt(n, n);
    vector<RealT, RealAlloc> s(std::min(m, n));

    char jobu = 'A';
    char jobvt = 'A';
    int ldu = n;
    int ldvt = m;
    int info = 0;

    T work_query;
    int lwork = -1;

    // Required rwork size for complex gesvd: 5 * min(m, n)
    constexpr int rwork_factor = 5;

    if constexpr (std::is_same_v<T, float>)
    {
      sgesvd_(&jobu, &jobvt, &lap_m, &lap_n, a_copy.begin(), &lda, s.begin(), vt.begin(), &ldu, u.begin(), &ldvt,
              &work_query, &lwork, &info);

      lwork = static_cast<int>(work_query);
      vector<T, Alloc> work(lwork);

      sgesvd_(&jobu, &jobvt, &lap_m, &lap_n, a_copy.begin(), &lda, s.begin(), vt.begin(), &ldu, u.begin(), &ldvt,
              work.begin(), &lwork, &info);
    }
    else if constexpr (std::is_same_v<T, double>)
    {
      dgesvd_(&jobu, &jobvt, &lap_m, &lap_n, a_copy.begin(), &lda, s.begin(), vt.begin(), &ldu, u.begin(), &ldvt,
              &work_query, &lwork, &info);

      lwork = static_cast<int>(work_query);
      vector<T, Alloc> work(lwork);

      dgesvd_(&jobu, &jobvt, &lap_m, &lap_n, a_copy.begin(), &lda, s.begin(), vt.begin(), &ldu, u.begin(), &ldvt,
              work.begin(), &lwork, &info);
    }
    else if constexpr (std::is_same_v<T, std::complex<float>>)
    {
      vector<RealT, RealAlloc> rwork(rwork_factor * std::min(m, n));

      cgesvd_(&jobu, &jobvt, &lap_m, &lap_n, a_copy.begin(), &lda, s.begin(), vt.begin(), &ldu, u.begin(), &ldvt,
              &work_query, &lwork, rwork.begin(), &info);

      lwork = static_cast<int>(work_query.real());
      vector<T, Alloc> work(lwork);

      cgesvd_(&jobu, &jobvt, &lap_m, &lap_n, a_copy.begin(), &lda, s.begin(), vt.begin(), &ldu, u.begin(), &ldvt,
              work.begin(), &lwork, rwork.begin(), &info);
    }
    else if constexpr (std::is_same_v<T, std::complex<double>>)
    {
      vector<RealT, RealAlloc> rwork(rwork_factor * std::min(m, n));

      zgesvd_(&jobu, &jobvt, &lap_m, &lap_n, a_copy.begin(), &lda, s.begin(), vt.begin(), &ldu, u.begin(), &ldvt,
              &work_query, &lwork, rwork.begin(), &info);

      lwork = static_cast<int>(work_query.real());
      vector<T, Alloc> work(lwork);

      zgesvd_(&jobu, &jobvt, &lap_m, &lap_n, a_copy.begin(), &lda, s.begin(), vt.begin(), &ldu, u.begin(), &ldvt,
              work.begin(), &lwork, rwork.begin(), &info);
    }

    if (info == 0)
    {
      return svd_result<T, Alloc>(std::move(u), std::move(s), std::move(vt));
    }
  }
#endif
  return std::nullopt;
}

export template<typename M>
  requires lam::concepts::experimental::matrix_c_weak<M, typename M::scalar_type>
auto pinv(const M& A) -> matrix<typename M::scalar_type> // NOLINT(readability-function-cognitive-complexity)
{
  using T = typename M::scalar_type;
  using Alloc = std::allocator<T>;
  using RealT = std::conditional_t<std::is_same_v<T, std::complex<double>>, double,
                                   std::conditional_t<std::is_same_v<T, std::complex<float>>, float, T>>;
  RealT tol = std::numeric_limits<RealT>::epsilon();

  auto svd_opt = svd(A);
  if (!svd_opt.has_value())
  {
    return matrix<T, Alloc>(A.cols(), A.rows());
  }

  auto& res = svd_opt.value();
  std::size_t m = A.rows();
  std::size_t n = A.cols();
  const std::size_t min_mn = std::min(m, n);

  RealT max_s = 0;
  for (std::size_t i = 0; i < min_mn; ++i)
  {
    if (res.s[i] > max_s)
    {
      max_s = res.s[i];
    }
  }

  RealT threshold = std::max(m, n) * max_s * std::abs(tol);

  matrix<T, Alloc> s_inv(n, m);
  for (std::size_t i = 0; i < n; ++i)
  {
    for (std::size_t j = 0; j < m; ++j)
    {
      s_inv[i, j] = T{0};
    }
  }

  for (std::size_t i = 0; i < min_mn; ++i)
  {
    if (res.s[i] > threshold)
    {
      s_inv[i, i] = T{1} / static_cast<T>(res.s[i]);
    }
  }

  // A = U S V^T  ==>  A^+ = V S^+ U^T
  // But wait! U is m x m, S is m x n, V^T is n x n.
  // Wait, S is formally an m x n diagonal matrix. And S^+ is n x m.
  // We have U (m x m), vt (n x n). So V is vt^T (which is n x n for complex, vt^H).

  matrix<T, Alloc> u_h(m, m);
  for (std::size_t i = 0; i < m; ++i)
  {
    for (std::size_t j = 0; j < m; ++j)
    {
      if constexpr (std::is_same_v<T, std::complex<double>> || std::is_same_v<T, std::complex<float>>)
      {
        u_h[i, j] = std::conj(res.u[j, i]);
      }
      else
      {
        u_h[i, j] = res.u[j, i];
      }
    }
  }

  matrix<T, Alloc> v(n, n);
  for (std::size_t i = 0; i < n; ++i)
  {
    for (std::size_t j = 0; j < n; ++j)
    {
      if constexpr (std::is_same_v<T, std::complex<double>> || std::is_same_v<T, std::complex<float>>)
      {
        v[i, j] = std::conj(res.vt[j, i]);
      }
      else
      {
        v[i, j] = res.vt[j, i];
      }
    }
  }

  auto temp = v * s_inv;
  return temp * u_h;
}

export template<typename M, typename V>
  requires lam::concepts::experimental::matrix_c_weak<M, typename M::scalar_type> &&
           lam::concepts::experimental::vector_c_weak<V>
auto lstsq(const M& A, const V& b) -> vector<typename M::scalar_type>
{
  using T = typename M::scalar_type;
  using Alloc = std::allocator<T>;

  matrix<T, Alloc> pinv_A = pinv(A);
  vector<T, Alloc> x(pinv_A.rows());

  for (std::size_t i = 0; i < pinv_A.rows(); ++i)
  {
    T sum = T{0};
    for (std::size_t j = 0; j < pinv_A.cols(); ++j)
    {
      sum += pinv_A[i, j] * b[j];
    }
    x[i] = sum;
  }

  return x;
}

export template<typename M>
  requires lam::concepts::experimental::matrix_c_weak<M, typename M::scalar_type>
std::size_t matrix_rank(const M& A)
{
  using T = typename M::scalar_type;
  using RealT = std::conditional_t<std::is_same_v<T, std::complex<double>>, double,
                                   std::conditional_t<std::is_same_v<T, std::complex<float>>, float, T>>;
  RealT tol = std::numeric_limits<RealT>::epsilon();

  auto svd_opt = svd(A);
  if (!svd_opt.has_value())
  {
    return 0;
  }

  auto& res = svd_opt.value();
  std::size_t m = A.rows();
  std::size_t n = A.cols();
  const std::size_t min_mn = std::min(m, n);

  RealT max_s = 0;
  for (std::size_t i = 0; i < min_mn; ++i)
  {
    if (res.s[i] > max_s)
    {
      max_s = res.s[i];
    }
  }

  RealT threshold = std::max(m, n) * max_s * std::abs(tol);

  std::size_t r = 0;
  for (std::size_t i = 0; i < min_mn; ++i)
  {
    if (res.s[i] > threshold)
    {
      r++;
    }
  }
  return r;
}

// Real scalar type underlying T: T for real, float/double for complex
template<typename T>
struct real_of
{
  using type = T;
};
template<typename T>
struct real_of<std::complex<T>>
{
  using type = T;
};
template<typename T>
using real_type_t = typename real_of<T>::type;

// Conjugate: identity for real types, std::conj for complex
// Non-exported module symbols already have module linkage — NOLINT(misc-use-internal-linkage)
template<typename T>
T conjugate_of(T x)
{ return x; } // NOLINT(misc-use-internal-linkage)
template<typename T>
std::complex<T> conjugate_of(std::complex<T> x)
{ return std::conj(x); } // NOLINT(misc-use-internal-linkage)

/**
 * Per-call SVD error diagnostics.
 *
 * reconstruction_error   — ||A − U·diag(s)·Vt||_F / ||A||_F  (backward error; O(ε) for a correct factorisation)
 * u_orthogonality_error  — ||U^H·U − I||_F                    (O(√m·ε))
 * vt_orthogonality_error — ||Vt·Vt^H − I||_F                  (O(√n·ε))
 * condition_number       — s[0]/s[k−1], or ∞ if rank-deficient
 * singular_value_bound   — ε·||A||_F  (Weyl's theorem: absolute error per singular value)
 * min_singular_gap       — min_i(s[i]−s[i+1]); singular vector perturbation scales as ε·||A||/gap
 *
 * Note: orthogonality checks are O(m³) and O(n³); intended for diagnostic/testing use.
 */
export template<typename T>
struct svd_diagnostics
{
  using real_type = real_type_t<T>;
  real_type reconstruction_error;
  real_type u_orthogonality_error;
  real_type vt_orthogonality_error;
  real_type condition_number;
  real_type singular_value_bound;
  real_type min_singular_gap;
};

/**
 * Compute SVD error diagnostics for the factorisation res of A = U·diag(s)·Vt.
 */
export template<typename M>
  requires lam::concepts::experimental::matrix_c_weak<M, typename M::scalar_type>
auto svd_check(const M& A,
               const svd_result<typename M::scalar_type>& res) // NOLINT(readability-function-cognitive-complexity)
  -> svd_diagnostics<typename M::scalar_type>
{
  using T = typename M::scalar_type;
  using RealT = real_type_t<T>;

  const std::size_t m = A.rows();
  const std::size_t n = A.cols();
  const std::size_t k = res.s.size();

  // Three independent branches run concurrently:
  //   - reconstruction error + ||A||_F²  (O(mn·k))
  //   - ||U^H·U − I||_F                  (O(m³))
  //   - ||Vt·Vt^H − I||_F                (O(n³))
  // All three only read from A and res; each writes only to its own local.

  auto fut_recon = std::async(std::launch::async, [&]() -> std::pair<RealT, RealT> {
    RealT frob_A_sq = RealT{0};
    RealT frob_res_sq = RealT{0};
    for (std::size_t i = 0; i < m; ++i)
    {
      for (std::size_t j = 0; j < n; ++j)
      {
        frob_A_sq += std::norm(A[i, j]);
        T rec = T{0};
        for (std::size_t r = 0; r < k; ++r)
        {
          rec += res.u[i, r] * static_cast<T>(res.s[r]) * res.vt[r, j];
        }
        frob_res_sq += std::norm(A[i, j] - rec);
      }
    }
    return {frob_A_sq, frob_res_sq};
  });

  auto fut_u = std::async(std::launch::async, [&]() -> RealT {
    RealT u_orth_sq = RealT{0};
    for (std::size_t i = 0; i < m; ++i)
    {
      for (std::size_t j = 0; j < m; ++j)
      {
        T dot = T{0};
        for (std::size_t r = 0; r < m; ++r)
        {
          dot += conjugate_of(res.u[r, i]) * res.u[r, j];
        }
        T expected = (i == j) ? T{1} : T{0};
        u_orth_sq += std::norm(dot - expected);
      }
    }
    return u_orth_sq;
  });

  auto fut_vt = std::async(std::launch::async, [&]() -> RealT {
    RealT vt_orth_sq = RealT{0};
    for (std::size_t i = 0; i < n; ++i)
    {
      for (std::size_t j = 0; j < n; ++j)
      {
        T dot = T{0};
        for (std::size_t r = 0; r < n; ++r)
        {
          dot += res.vt[i, r] * conjugate_of(res.vt[j, r]);
        }
        T expected = (i == j) ? T{1} : T{0};
        vt_orth_sq += std::norm(dot - expected);
      }
    }
    return vt_orth_sq;
  });

  auto [frob_A_sq, frob_res_sq] = fut_recon.get();
  RealT u_orth_sq = fut_u.get();
  RealT vt_orth_sq = fut_vt.get();

  RealT reconstruction_error = (frob_A_sq > RealT{0}) ? std::sqrt(frob_res_sq / frob_A_sq) : std::sqrt(frob_res_sq);

  // Condition number
  RealT cond = (k > 0 && res.s[k - 1] > RealT{0}) ? res.s[0] / res.s[k - 1] : std::numeric_limits<RealT>::infinity();

  // Weyl bound: absolute error per singular value ≤ ε·||A||_F
  RealT sv_bound = std::sqrt(frob_A_sq) * std::numeric_limits<RealT>::epsilon();

  // Minimum consecutive singular value gap (bounds singular vector sensitivity)
  RealT min_gap = std::numeric_limits<RealT>::infinity();
  for (std::size_t i = 0; i + 1 < k; ++i)
  {
    min_gap = std::min(min_gap, res.s[i] - res.s[i + 1]);
  }

  return {reconstruction_error, std::sqrt(u_orth_sq), std::sqrt(vt_orth_sq), cond, sv_bound, min_gap};
}

} // namespace lam::linalg
