/*
 *  matrix_eigenvalue.cppm - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  matrix_eigenvalue is a c++ module providing eigenvalue
 *  decomposition using LAPACK via Accelerate/OpenBLAS
 */

export module lam.linearalgebra:matrix.eigenvalue;

import std;
import :config;
import :vectorspace;
import :matrix;

namespace lam::linalg
{

template<class>
inline constexpr bool always_false_v = false;

namespace detail
{
// ||M||_F. Every tolerance here is measured against it rather than assuming entries of order one.
template<typename M>
double frobenius_norm(const M& m)
{
  double sum = 0.0;
  for (std::size_t i = 0; i < m.rows(); ++i)
    for (std::size_t j = 0; j < m.cols(); ++j)
      sum += std::norm(m[i, j]);
  return std::sqrt(sum);
}
} // namespace detail

/**
 * What a generalized eigenvalue actually is.
 *
 * A generalized eigenvalue is the ratio alpha/beta returned by the pencil solver. A negligible
 * beta means the ratio is not finite; a negligible alpha as well means it is not determined at
 * all (0/0), and the reported value carries no information. Callers deciding on an eigenvalue
 * must consult this rather than re-testing the value.
 */
export enum class eigenvalue_kind : std::uint8_t
{
  finite,
  infinite,
  indeterminate
};

#ifdef LAM_USE_BLAS
extern "C"
{
  void dsygv_(const int* itype, const char* jobz, const char* uplo, const int* n, double* a, const int* lda, double* b,
              const int* ldb, double* w, double* work, const int* lwork, int* info);

  void dsygvd_(const int* itype, const char* jobz, const char* uplo, const int* n, double* a, const int* lda, double* b,
               const int* ldb, double* w, double* work, const int* lwork, int* iwork, const int* liwork, int* info);

  void dsyev_(const char* jobz, const char* uplo, const int* n, double* a, const int* lda, double* w, double* work,
              const int* lwork, int* info);

  void zggev_(const char* jobvl, const char* jobvr, const int* n, std::complex<double>* a, const int* lda,
              std::complex<double>* b, const int* ldb, std::complex<double>* alpha, std::complex<double>* beta,
              std::complex<double>* vl, const int* ldvl, std::complex<double>* vr, const int* ldvr,
              std::complex<double>* work, const int* lwork, double* rwork, int* info);

  void zhegvd_(const int* itype, const char* jobz, const char* uplo, const int* n, std::complex<double>* a,
               const int* lda, std::complex<double>* b, const int* ldb, double* w, std::complex<double>* work,
               const int* lwork, double* rwork, const int* lrwork, int* iwork, const int* liwork, int* info);

  void zhegv_(const int* itype, const char* jobz, const char* uplo, const int* n, std::complex<double>* a,
              const int* lda, std::complex<double>* b, const int* ldb, double* w, std::complex<double>* work,
              const int* lwork, double* rwork, int* info);

  void zgeev_(const char* jobvl, const char* jobvr, const int* n, std::complex<double>* a, const int* lda,
              std::complex<double>* w, std::complex<double>* vl, const int* ldvl, std::complex<double>* vr,
              const int* ldvr, std::complex<double>* work, const int* lwork, double* rwork, int* info);
}
#endif

/**
 * Result of eigenvalue decomposition.
 */
export template<typename T, typename Alloc = std::allocator<T>>
struct eigen_result
{
  vector<T, Alloc> eigenvalues;                             // Eigenvalues in ascending order
  matrix<T, storage_layout::row_major, Alloc> eigenvectors; // Column k is eigenvector for eigenvalue k
  bool success;
};

/**
 * Result of generalized eigenvalue problem.
 * Solves: A * v = λ * B * v
 * where B must be symmetric positive definite.
 */
export template<typename T, typename Alloc = std::allocator<T>>
struct gep_result
{
  using kind_alloc = typename std::allocator_traits<Alloc>::template rebind_alloc<eigenvalue_kind>;

  vector<T, Alloc> eigenvalues;                             // Eigenvalues λ_k in ascending order
  matrix<T, storage_layout::row_major, Alloc> eigenvectors; // Column k is eigenvector for λ_k, v^T * B * v = 1
  bool success;
  // Parallel to eigenvalues; empty only when the solve failed before classification.
  std::vector<eigenvalue_kind, kind_alloc> kinds{};
};

/**
 * Solve generalized eigenvalue problem: A * v = λ * B * v
 *
 * Supports double (using dsygv) and std::complex<double> (using zggev).
 *
 * @param A n×n matrix
 * @param B n×n matrix (symmetric/hermitian positive definite for real case)
 * @return gep_result with eigenvalues and eigenvectors
 */
export template<typename M>
auto solve_gep(const M& A, const M& B) -> gep_result<typename M::scalar_type, detail::result_alloc_t<M>>
{
  using T = typename M::scalar_type;
  using Alloc = detail::result_alloc_t<M>;
  using RealAlloc = typename std::allocator_traits<Alloc>::template rebind_alloc<double>;
  using IntAlloc = typename std::allocator_traits<Alloc>::template rebind_alloc<int>;
  using kind_alloc = typename gep_result<T, Alloc>::kind_alloc;
  const Alloc alloc = detail::result_allocator(A);
  const std::size_t n = A.rows();

  // Both Hermitian routines require B positive definite, so every eigenvalue they return is finite.
  auto all_finite = [&] {
    return std::vector<eigenvalue_kind, kind_alloc>(n, eigenvalue_kind::finite, kind_alloc(alloc));
  };

  if (n == 0 || A.cols() != n || B.rows() != n || B.cols() != n)
  {
    return {vector<T, Alloc>(alloc), matrix<T, storage_layout::row_major, Alloc>(alloc), false};
  }

#ifdef LAM_USE_BLAS
  if constexpr (std::is_same_v<T, double>)
  {
    // Real Symmetric Definite GEP (dsygv/dsygvd)
    vector<T, Alloc> a_data(n * n, alloc);
    vector<T, Alloc> b_data(n * n, alloc);

    for (std::size_t j = 0; j < n; ++j)
    {
      for (std::size_t i = 0; i < n; ++i)
      {
        a_data[i + j * n] = A[i, j];
        b_data[i + j * n] = B[i, j];
      }
    }

    vector<T, Alloc> eigenvalues(n, alloc);
    int itype = 1;
    char jobz = 'V';
    char uplo = 'U';
    int n_int = static_cast<int>(n);
    int lda = n_int;
    int ldb = n_int;
    int info = 0;

    // Use Divide and Conquer (dsygvd) for larger matrices.
    // Crossover at n=100: dsygvd's O(n^2) memory and setup cost outweigh its
    // O(n^2) vs O(n^3) eigenvector advantage for small matrices. This threshold
    // is platform-dependent (cache hierarchy, LAPACK implementation).
    if (n > 100)
    {
      int lwork = -1;
      int liwork = -1;
      T work_query{};
      int iwork_query{};

      dsygvd_(&itype, &jobz, &uplo, &n_int, a_data.begin(), &lda, b_data.begin(), &ldb, eigenvalues.begin(),
              &work_query, &lwork, &iwork_query, &liwork, &info);

      const auto lw = detail::checked_lwork(work_query, info);
      const auto liw = detail::checked_lwork(iwork_query, info);
      if (!lw || !liw)
        return {vector<T, Alloc>(alloc), matrix<T, storage_layout::row_major, Alloc>(alloc), false};
      lwork = *lw;
      liwork = *liw;

      vector<T, Alloc> work(static_cast<std::size_t>(lwork), alloc);
      vector<int, IntAlloc> iwork(static_cast<std::size_t>(liwork), IntAlloc(alloc));

      dsygvd_(&itype, &jobz, &uplo, &n_int, a_data.begin(), &lda, b_data.begin(), &ldb, eigenvalues.begin(),
              work.begin(), &lwork, iwork.begin(), &liwork, &info);
    }
    else
    {
      int lwork = -1;
      T work_query{};
      dsygv_(&itype, &jobz, &uplo, &n_int, a_data.begin(), &lda, b_data.begin(), &ldb, eigenvalues.begin(), &work_query,
             &lwork, &info);

      const auto lw = detail::checked_lwork(work_query, info);
      if (!lw)
        return {vector<T, Alloc>(alloc), matrix<T, storage_layout::row_major, Alloc>(alloc), false};
      lwork = *lw;
      vector<T, Alloc> work(static_cast<std::size_t>(lwork), alloc);

      dsygv_(&itype, &jobz, &uplo, &n_int, a_data.begin(), &lda, b_data.begin(), &ldb, eigenvalues.begin(),
             work.begin(), &lwork, &info);
    }

    if (info != 0)
    {
      return {vector<T, Alloc>(alloc), matrix<T, storage_layout::row_major, Alloc>(alloc), false};
    }

    matrix<T, storage_layout::row_major, Alloc> eigenvectors(n, n, alloc);
    for (std::size_t j = 0; j < n; ++j)
    {
      for (std::size_t i = 0; i < n; ++i)
      {
        eigenvectors[i, j] = a_data[i + j * n];
      }
    }
    return {std::move(eigenvalues), std::move(eigenvectors), true, all_finite()};
  }
  else if constexpr (std::is_same_v<T, std::complex<double>>)
  {
    vector<T, Alloc> a_data(n * n, alloc);
    vector<T, Alloc> b_data(n * n, alloc);

    for (std::size_t j = 0; j < n; ++j)
    {
      for (std::size_t i = 0; i < n; ++i)
      {
        a_data[i + j * n] = A[i, j];
        b_data[i + j * n] = B[i, j];
      }
    }

    vector<T, Alloc> eigenvalues(n, alloc);
    int n_int = static_cast<int>(n);
    int lda = n_int;
    int ldb = n_int;

    // Check if A is Hermitian before attempting Hermitian solver —
    // zhegvd_ silently reads only the upper triangle and produces wrong
    // results for non-Hermitian input (without setting info != 0).
    // Conservative: a wrong "yes" feeds zhegvd_ garbage, a wrong "no" only costs the general routine.
    const double hermitian_tol =
      static_cast<double>(n) * std::numeric_limits<double>::epsilon() * detail::frobenius_norm(A);
    bool a_hermitian = true;
    for (std::size_t i = 0; i < n && a_hermitian; ++i)
    {
      // A Hermitian matrix has a real diagonal; zhegvd_ assumes it without checking.
      if (std::abs(A[i, i].imag()) > hermitian_tol)
        a_hermitian = false;
      for (std::size_t j = i + 1; j < n && a_hermitian; ++j)
        if (std::abs(A[i, j] - std::conj(A[j, i])) > hermitian_tol)
          a_hermitian = false;
    }

    int hermitian_info = -1; // -1 = not attempted
    if (a_hermitian)
    {
      int itype = 1;
      char jobz = 'V';
      char uplo = 'U';

      int lwork = -1;
      int lrwork = -1;
      int liwork = -1;
      std::complex<double> work_query{};
      double rwork_query{};
      int iwork_query{};

      zhegvd_(&itype, &jobz, &uplo, &n_int, a_data.begin(), &lda, b_data.begin(), &ldb,
              reinterpret_cast<double*>(eigenvalues.begin()), &work_query, &lwork, &rwork_query, &lrwork, &iwork_query,
              &liwork, &hermitian_info);

      const auto lw = detail::checked_lwork(work_query.real(), hermitian_info);
      const auto lrw = detail::checked_lwork(rwork_query, hermitian_info);
      const auto liw = detail::checked_lwork(iwork_query, hermitian_info);
      // A failed query counts as zhegvd_ failing, so the zggev_ fallback below takes over.
      if (!lw || !lrw || !liw)
        hermitian_info = -1;
      else
      {
        lwork = *lw;
        lrwork = *lrw;
        liwork = *liw;

        vector<std::complex<double>, Alloc> work(static_cast<std::size_t>(lwork), alloc);
        vector<double, RealAlloc> rwork(static_cast<std::size_t>(lrwork), RealAlloc(alloc));
        vector<int, IntAlloc> iwork(static_cast<std::size_t>(liwork), IntAlloc(alloc));

        vector<double, RealAlloc> real_eigenvalues(n, RealAlloc(alloc));

        zhegvd_(&itype, &jobz, &uplo, &n_int, a_data.begin(), &lda, b_data.begin(), &ldb, real_eigenvalues.begin(),
                work.begin(), &lwork, rwork.begin(), &lrwork, iwork.begin(), &liwork, &hermitian_info);

        if (hermitian_info == 0)
        {
          for (std::size_t i = 0; i < n; ++i)
            eigenvalues[i] = real_eigenvalues[i];
        }
      }
    }

    // Fallback to zggev if A is not Hermitian or zhegvd_ failed
    if (hermitian_info != 0)
    {
      // Re-copy A and B (zhegvd_ destroys them on failure)
      if constexpr (M::layout == storage_layout::col_major)
      {
        std::copy(A.begin(), A.end(), a_data.begin());
        std::copy(B.begin(), B.end(), b_data.begin());
      }
      else
      {
        for (std::size_t j = 0; j < n; ++j)
        {
          for (std::size_t i = 0; i < n; ++i)
          {
            a_data[i + j * n] = A[i, j];
            b_data[i + j * n] = B[i, j];
          }
        }
      }

      int info = 0;
      vector<T, Alloc> alpha(n, alloc);
      vector<T, Alloc> beta(n, alloc);
      vector<T, Alloc> vl(n * n, alloc);
      vector<T, Alloc> vr(n * n, alloc);
      vector<double, RealAlloc> rwork(8 * n, RealAlloc(alloc));

      char jobvl = 'N';
      char jobvr = 'V';
      int ldvl = n_int;
      int ldvr = n_int;

      T work_query{};
      int lwork = -1;
      zggev_(&jobvl, &jobvr, &n_int, a_data.begin(), &lda, b_data.begin(), &ldb, alpha.begin(), beta.begin(),
             vl.begin(), &ldvl, vr.begin(), &ldvr, &work_query, &lwork, rwork.begin(), &info);

      const auto lw = detail::checked_lwork(work_query.real(), info);
      if (!lw)
        return {vector<T, Alloc>(alloc), matrix<T, storage_layout::row_major, Alloc>(alloc), false};
      lwork = *lw;
      vector<T, Alloc> work(static_cast<std::size_t>(lwork), alloc);

      zggev_(&jobvl, &jobvr, &n_int, a_data.begin(), &lda, b_data.begin(), &ldb, alpha.begin(), beta.begin(),
             vl.begin(), &ldvl, vr.begin(), &ldvr, work.begin(), &lwork, rwork.begin(), &info);

      if (info != 0)
      {
        return {vector<T, Alloc>(alloc), matrix<T, storage_layout::row_major, Alloc>(alloc), false};
      }

      struct eigen_pair
      {
        T eigenvalue;
        std::size_t index;
        eigenvalue_kind kind;
      };

      std::vector<eigen_pair, typename std::allocator_traits<Alloc>::template rebind_alloc<eigen_pair>> pairs(
        n, typename std::allocator_traits<Alloc>::template rebind_alloc<eigen_pair>(alloc));

      // alpha and beta come out of zggev_ scaled to the pencil, so the thresholds are too.
      const double beta_tol =
        static_cast<double>(n) * std::numeric_limits<double>::epsilon() * detail::frobenius_norm(B);
      const double alpha_tol =
        static_cast<double>(n) * std::numeric_limits<double>::epsilon() * detail::frobenius_norm(A);

      for (std::size_t i = 0; i < n; ++i)
      {
        const bool beta_negligible = std::abs(beta[i]) <= beta_tol;
        if (beta_negligible && std::abs(alpha[i]) <= alpha_tol)
          pairs[i] = {T{std::numeric_limits<double>::quiet_NaN(), 0.0}, i, eigenvalue_kind::indeterminate};
        else if (beta_negligible)
          pairs[i] = {T{std::numeric_limits<double>::infinity(), 0.0}, i, eigenvalue_kind::infinite};
        else
          pairs[i] = {alpha[i] / beta[i], i, eigenvalue_kind::finite};
      }

      // Kind first: an indeterminate value is NaN, and comparing it is not a strict weak ordering.
      std::sort(pairs.begin(), pairs.end(), [](const eigen_pair& a, const eigen_pair& b) {
        if (a.kind != b.kind)
          return a.kind < b.kind;
        if (a.kind != eigenvalue_kind::finite)
          return false;
        return a.eigenvalue.real() < b.eigenvalue.real();
      });

      std::vector<eigenvalue_kind, kind_alloc> kinds(n, eigenvalue_kind::finite, kind_alloc(alloc));

      matrix<T, storage_layout::row_major, Alloc> eigenvectors(n, n, alloc);
      for (std::size_t i = 0; i < n; ++i)
      {
        eigenvalues[i] = pairs[i].eigenvalue;
        kinds[i] = pairs[i].kind;
        std::size_t orig_idx = pairs[i].index;

        for (std::size_t r = 0; r < n; ++r)
          eigenvectors[r, i] = vr[r + orig_idx * n];

        // Enforce B-normalization: v^H * B * v = 1
        T norm_sq = T(0);
        for (std::size_t r = 0; r < n; ++r)
        {
          T B_row_dot_v = T(0);
          for (std::size_t c = 0; c < n; ++c)
          {
            B_row_dot_v += B[r, c] * eigenvectors[c, i];
          }
          norm_sq += std::conj(eigenvectors[r, i]) * B_row_dot_v;
        }

        // v^H B v is a quadratic form in B, so its scale is ||B||, not one.
        if (std::abs(norm_sq) > beta_tol * detail::frobenius_norm(B))
        {
          double scale = 1.0 / std::sqrt(std::abs(norm_sq));
          for (std::size_t r = 0; r < n; ++r)
          {
            eigenvectors[r, i] *= scale;
          }
        }
      }
      return {std::move(eigenvalues), std::move(eigenvectors), true, std::move(kinds)};
    }

    matrix<T, storage_layout::row_major, Alloc> eigenvectors(n, n, alloc);
    for (std::size_t j = 0; j < n; ++j)
    {
      for (std::size_t i = 0; i < n; ++i)
      {
        eigenvectors[i, j] = a_data[i + j * n];
      }
    }
    return {std::move(eigenvalues), std::move(eigenvectors), true, all_finite()};
  }
  else
  {
    static_assert(always_false_v<T>, "solve_gep not implemented for this scalar type");
    return {vector<T, Alloc>(alloc), matrix<T, storage_layout::row_major, Alloc>(alloc), false};
  }
#else
  return {vector<T, Alloc>(alloc), matrix<T, storage_layout::row_major, Alloc>(alloc), false};
#endif
}

/**
 * Compute eigenvalues and eigenvectors of a symmetric matrix.
 *
 * @param A n×n symmetric matrix (only upper triangle is read)
 * @return eigen_result with eigenvalues and eigenvectors
 */
export template<typename M>
  requires std::same_as<typename M::scalar_type, double>
auto symmetric_eigen(const M& A) -> eigen_result<double, detail::result_alloc_t<M>>
{
  using T = double;
  using Alloc = detail::result_alloc_t<M>;
  using RealAlloc = typename std::allocator_traits<Alloc>::template rebind_alloc<double>;
  using IntAlloc = typename std::allocator_traits<Alloc>::template rebind_alloc<int>;
  const Alloc alloc = detail::result_allocator(A);

  const std::size_t n = A.rows();

  if (n == 0 || A.cols() != n)
  {
    return {vector<T, Alloc>(alloc), matrix<T, storage_layout::row_major, Alloc>(alloc), false};
  }

#ifdef LAM_USE_BLAS
  // Create column-major copy for LAPACK
  vector<T, Alloc> a_data(n * n, alloc);

  if constexpr (M::layout == storage_layout::col_major)
  {
    std::copy(A.begin(), A.end(), a_data.begin());
  }
  else
  {
    for (std::size_t j = 0; j < n; ++j)
    {
      for (std::size_t i = 0; i < n; ++i)
      {
        a_data[i + j * n] = A[i, j];
      }
    }
  }

  vector<T, Alloc> eigenvalues(n, alloc);
  char jobz = 'V';
  char uplo = 'U';
  int n_int = static_cast<int>(n);
  int lda = n_int;
  int info = 0;

  // Query optimal workspace size
  int lwork = -1;
  T work_query{};
  dsyev_(&jobz, &uplo, &n_int, a_data.begin(), &lda, eigenvalues.begin(), &work_query, &lwork, &info);

  const auto lw = detail::checked_lwork(work_query, info);
  if (!lw)
    return {vector<T, Alloc>(alloc), matrix<T, storage_layout::row_major, Alloc>(alloc), false};
  lwork = *lw;
  vector<T, Alloc> work(static_cast<std::size_t>(lwork), alloc);

  // Compute eigenvalues and eigenvectors
  dsyev_(&jobz, &uplo, &n_int, a_data.begin(), &lda, eigenvalues.begin(), work.begin(), &lwork, &info);

  if (info != 0)
  {
    return {vector<T, Alloc>(alloc), matrix<T, storage_layout::row_major, Alloc>(alloc), false};
  }

  // Copy eigenvectors to matrix
  matrix<T, storage_layout::row_major, Alloc> eigenvectors(n, n, alloc);
  for (std::size_t j = 0; j < n; ++j)
  {
    for (std::size_t i = 0; i < n; ++i)
    {
      eigenvectors[i, j] = a_data[i + j * n];
    }
  }

  return {std::move(eigenvalues), std::move(eigenvectors), true};
#else
  return {vector<T, Alloc>(alloc), matrix<T, storage_layout::row_major, Alloc>(alloc), false};
#endif
}

/**
 * Compute eigenvalues and eigenvectors of a general complex matrix.
 *
 * Solves the standard eigenvalue problem A * v = λ * v for arbitrary
 * (non-Hermitian) complex matrices using LAPACK's zgeev.
 *
 * @param A n×n complex matrix
 * @return eigen_result with complex eigenvalues and eigenvectors
 */
export template<typename M>
  requires std::same_as<typename M::scalar_type, std::complex<double>>
auto general_complex_eigen(const M& A) -> eigen_result<std::complex<double>, detail::result_alloc_t<M>>
{
  using T = std::complex<double>;
  using Alloc = detail::result_alloc_t<M>;
  using RealAlloc = typename std::allocator_traits<Alloc>::template rebind_alloc<double>;
  using IntAlloc = typename std::allocator_traits<Alloc>::template rebind_alloc<int>;
  const Alloc alloc = detail::result_allocator(A);

  const std::size_t n = A.rows();

  if (n == 0 || A.cols() != n)
  {
    return {vector<T, Alloc>(alloc), matrix<T, storage_layout::row_major, Alloc>(alloc), false};
  }

#ifdef LAM_USE_BLAS
  vector<T, Alloc> a_data(n * n, alloc);

  if constexpr (M::layout == storage_layout::col_major)
  {
    std::copy(A.begin(), A.end(), a_data.begin());
  }
  else
  {
    for (std::size_t j = 0; j < n; ++j)
    {
      for (std::size_t i = 0; i < n; ++i)
      {
        a_data[i + j * n] = A[i, j];
      }
    }
  }

  vector<T, Alloc> eigenvalues(n, alloc);
  int n_int = static_cast<int>(n);
  int lda = n_int;
  int info = 0;

  char jobvl = 'N';
  char jobvr = 'V';
  int ldvl = 1;
  int ldvr = n_int;

  vector<T, Alloc> vl(1, alloc);
  vector<T, Alloc> vr(n * n, alloc);
  vector<double, RealAlloc> rwork(2 * n, RealAlloc(alloc));

  // Workspace query
  int lwork = -1;
  T work_query{};
  zgeev_(&jobvl, &jobvr, &n_int, a_data.begin(), &lda, eigenvalues.begin(), vl.begin(), &ldvl, vr.begin(), &ldvr,
         &work_query, &lwork, rwork.begin(), &info);

  const auto lw = detail::checked_lwork(work_query.real(), info);
  if (!lw)
    return {vector<T, Alloc>(alloc), matrix<T, storage_layout::row_major, Alloc>(alloc), false};
  lwork = *lw;
  vector<T, Alloc> work(static_cast<std::size_t>(lwork), alloc);

  zgeev_(&jobvl, &jobvr, &n_int, a_data.begin(), &lda, eigenvalues.begin(), vl.begin(), &ldvl, vr.begin(), &ldvr,
         work.begin(), &lwork, rwork.begin(), &info);

  if (info != 0)
  {
    return {vector<T, Alloc>(alloc), matrix<T, storage_layout::row_major, Alloc>(alloc), false};
  }

  matrix<T, storage_layout::row_major, Alloc> eigenvectors(n, n, alloc);
  for (std::size_t j = 0; j < n; ++j)
  {
    for (std::size_t i = 0; i < n; ++i)
    {
      eigenvectors[i, j] = vr[i + j * n];
    }
  }

  return {std::move(eigenvalues), std::move(eigenvectors), true};
#else
  return {vector<T, Alloc>(alloc), matrix<T, storage_layout::row_major, Alloc>(alloc), false};
#endif
}

} // namespace lam::linalg
