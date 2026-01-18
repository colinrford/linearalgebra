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
}
#endif

/**
 * Result of eigenvalue decomposition.
 */
export template<typename T>
struct eigen_result
{
  vector<T> eigenvalues;  // Eigenvalues in ascending order
  matrix<T> eigenvectors; // Column k is eigenvector for eigenvalue k
  bool success;
};

/**
 * Result of generalized eigenvalue problem.
 * Solves: A * v = λ * B * v
 * where B must be symmetric positive definite.
 */
export template<typename T>
struct gep_result
{
  vector<T> eigenvalues;  // Eigenvalues λ_k in ascending order
  matrix<T> eigenvectors; // Column k is eigenvector for λ_k, normalized so v^T * B * v = 1
  bool success;
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
auto solve_gep(const M& A, const M& B) -> gep_result<typename M::scalar_type>
{
  using T = typename M::scalar_type;
  const std::size_t n = A.rows();

  if (n == 0 || A.cols() != n || B.rows() != n || B.cols() != n)
  {
    return {vector<T>{}, matrix<T>{}, false};
  }

#ifdef LAM_USE_BLAS
  if constexpr (std::is_same_v<T, double>)
  {
    // Real Symmetric Definite GEP (dsygv/dsygvd)
    vector<T> a_data(n * n);
    vector<T> b_data(n * n);

    for (std::size_t j = 0; j < n; ++j)
    {
      for (std::size_t i = 0; i < n; ++i)
      {
        a_data[i + j * n] = A[i, j];
        b_data[i + j * n] = B[i, j];
      }
    }

    vector<T> eigenvalues(n);
    int itype = 1;
    char jobz = 'V';
    char uplo = 'U';
    int n_int = static_cast<int>(n);
    int lda = n_int;
    int ldb = n_int;
    int info = 0;

    // Use Divide and Conquer (dsygvd) for larger matrices (n > 100)
    if (n > 100)
    {
      int lwork = -1;
      int liwork = -1;
      T work_query;
      int iwork_query;

      dsygvd_(&itype, &jobz, &uplo, &n_int, a_data.begin(), &lda, b_data.begin(), &ldb, eigenvalues.begin(),
              &work_query, &lwork, &iwork_query, &liwork, &info);

      lwork = static_cast<int>(work_query);
      liwork = iwork_query;

      vector<T> work(static_cast<std::size_t>(lwork));
      vector<int> iwork(static_cast<std::size_t>(liwork));

      dsygvd_(&itype, &jobz, &uplo, &n_int, a_data.begin(), &lda, b_data.begin(), &ldb, eigenvalues.begin(),
              work.begin(), &lwork, iwork.begin(), &liwork, &info);
    }
    else
    {
      int lwork = -1;
      T work_query;
      dsygv_(&itype, &jobz, &uplo, &n_int, a_data.begin(), &lda, b_data.begin(), &ldb, eigenvalues.begin(), &work_query,
             &lwork, &info);

      lwork = static_cast<int>(work_query);
      vector<T> work(static_cast<std::size_t>(lwork));

      dsygv_(&itype, &jobz, &uplo, &n_int, a_data.begin(), &lda, b_data.begin(), &ldb, eigenvalues.begin(),
             work.begin(), &lwork, &info);
    }

    if (info != 0)
    {
      return {vector<T>{}, matrix<T>{}, false};
    }

    matrix<T> eigenvectors(n, n);
    for (std::size_t j = 0; j < n; ++j)
    {
      for (std::size_t i = 0; i < n; ++i)
      {
        eigenvectors[i, j] = a_data[i + j * n];
      }
    }
    return {std::move(eigenvalues), std::move(eigenvectors), true};
  }
  else if constexpr (std::is_same_v<T, std::complex<double>>)
  { // Complex General GEP (zggev)
    // Note: A and B are copied to column-major format
    vector<T> a_data(n * n);
    vector<T> b_data(n * n);

    for (std::size_t j = 0; j < n; ++j)
    {
      for (std::size_t i = 0; i < n; ++i)
      {
        a_data[i + j * n] = A[i, j];
        b_data[i + j * n] = B[i, j];
      }
    }

    vector<T> eigenvalues(n);
    int itype = 1;
    char jobz = 'V';
    char uplo = 'U';
    int n_int = static_cast<int>(n);
    int lda = n_int;
    int ldb = n_int;
    int info = 0;

    // Use Divide and Conquer (zhegvd) for larger matrices (n > 100)
    // Try Hermitian solver first
    int zhegvd_info = 0;
    if (n > 100)
    {
      int lwork = -1;
      std::complex<double> work_query;
      int lrwork = -1;
      double rwork_query;
      int liwork = -1;
      int iwork_query;

      zhegvd_(&itype, &jobz, &uplo, &n_int, a_data.begin(), &lda, b_data.begin(), &ldb,
              reinterpret_cast<double*>(eigenvalues.begin()), &work_query, &lwork, &rwork_query, &lrwork, &iwork_query,
              &liwork, &zhegvd_info);

      lwork = static_cast<int>(work_query.real());
      lrwork = static_cast<int>(rwork_query);
      liwork = iwork_query;

      vector<std::complex<double>> work(static_cast<std::size_t>(lwork));
      vector<double> rwork(static_cast<std::size_t>(lrwork));
      vector<int> iwork(static_cast<std::size_t>(liwork));

      vector<double> real_eigenvalues(n);

      zhegvd_(&itype, &jobz, &uplo, &n_int, a_data.begin(), &lda, b_data.begin(), &ldb, real_eigenvalues.begin(),
              work.begin(), &lwork, rwork.begin(), &lrwork, iwork.begin(), &liwork, &zhegvd_info);

      if (zhegvd_info == 0)
      {
        for (std::size_t i = 0; i < n; ++i)
          eigenvalues[i] = real_eigenvalues[i];
      }
    }
    else
    { // Force fallback for small n
      zhegvd_info = 1;
    }
    // Fallback to zggev (General GEP) if zhegvd failed (e.g. B not positive definite)
    // or if n <= 100 (where we didn't try zhegvd, or prefer QR/General)
    // Note: zhegvd destroys A and B on failure, so we must re-copy them.
    if (zhegvd_info != 0)
    { // Re-copy A and B
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

      vector<T> alpha(n);
      vector<T> beta(n);
      vector<T> vl(n * n); // Left (unused)
      vector<T> vr(n * n); // Right (used)
      vector<double> rwork(8 * n);

      char jobvl = 'N';
      char jobvr = 'V';
      int ldvl = n_int;
      int ldvr = n_int;

      T work_query;
      int lwork = -1;
      zggev_(&jobvl, &jobvr, &n_int, a_data.begin(), &lda, b_data.begin(), &ldb, alpha.begin(), beta.begin(),
             vl.begin(), &ldvl, vr.begin(), &ldvr, &work_query, &lwork, rwork.begin(), &info);

      lwork = static_cast<int>(work_query.real());
      vector<T> work(static_cast<std::size_t>(lwork));

      zggev_(&jobvl, &jobvr, &n_int, a_data.begin(), &lda, b_data.begin(), &ldb, alpha.begin(), beta.begin(),
             vl.begin(), &ldvl, vr.begin(), &ldvr, work.begin(), &lwork, rwork.begin(), &info);

      if (info != 0)
      {
        return {vector<T>{}, matrix<T>{}, false};
      }

      // Process generalized eigenvalues alpha/beta -> lambda
      // And sort them, because zggev returns them in random order
      struct eigen_pair
      {
        T eigenvalue;
        std::size_t index;
      };

      std::vector<eigen_pair> pairs(n);

      for (std::size_t i = 0; i < n; ++i)
      {
        if (std::abs(beta[i]) < 1e-14)
          pairs[i] = {std::numeric_limits<double>::infinity(), i};
        else
          pairs[i] = {alpha[i] / beta[i], i};
      }

      std::sort(pairs.begin(), pairs.end(),
                [](const eigen_pair& a, const eigen_pair& b) { return a.eigenvalue.real() < b.eigenvalue.real(); });

      matrix<T> eigenvectors(n, n);
      for (std::size_t i = 0; i < n; ++i)
      {
        eigenvalues[i] = pairs[i].eigenvalue;
        std::size_t orig_idx = pairs[i].index;

        for (std::size_t r = 0; r < n; ++r)
          eigenvectors[r, i] = vr[r + orig_idx * n];

        // Enforce B-normalization: v^H * B * v = 1
        // Compute dot = v^H * B * v
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

        if (std::abs(norm_sq) > 1e-14)
        {
          double scale = 1.0 / std::sqrt(std::abs(norm_sq));
          for (std::size_t r = 0; r < n; ++r)
          {
            eigenvectors[r, i] *= scale;
          }
        }
      }
      return {std::move(eigenvalues), std::move(eigenvectors), true};
    }

    matrix<T> eigenvectors(n, n);
    for (std::size_t j = 0; j < n; ++j)
    {
      for (std::size_t i = 0; i < n; ++i)
      {
        eigenvectors[i, j] = a_data[i + j * n];
      }
    }
    return {std::move(eigenvalues), std::move(eigenvectors), true};

    if (info != 0)
    {
      return {vector<T>{}, matrix<T>{}, false};
    }
  }
  else
  {
    static_assert(always_false_v<T>, "solve_gep not implemented for this scalar type");
    return {vector<T>{}, matrix<T>{}, false};
  }
#else
  return {vector<T>{}, matrix<T>{}, false};
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
auto symmetric_eigen(const M& A) -> eigen_result<double>
{
  using T = double;

  const std::size_t n = A.rows();

  if (n == 0 || A.cols() != n)
  {
    return {vector<T>{}, matrix<T>{}, false};
  }

#ifdef LAM_USE_BLAS
  // Create column-major copy for LAPACK
  vector<T> a_data(n * n);

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

  vector<T> eigenvalues(n);
  char jobz = 'V';
  char uplo = 'U';
  int n_int = static_cast<int>(n);
  int lda = n_int;
  int info = 0;

  // Query optimal workspace size
  int lwork = -1;
  T work_query;
  dsyev_(&jobz, &uplo, &n_int, a_data.begin(), &lda, eigenvalues.begin(), &work_query, &lwork, &info);

  lwork = static_cast<int>(work_query);
  vector<T> work(static_cast<std::size_t>(lwork));

  // Compute eigenvalues and eigenvectors
  dsyev_(&jobz, &uplo, &n_int, a_data.begin(), &lda, eigenvalues.begin(), work.begin(), &lwork, &info);

  if (info != 0)
  {
    return {vector<T>{}, matrix<T>{}, false};
  }

  // Copy eigenvectors to matrix
  matrix<T> eigenvectors(n, n);
  for (std::size_t j = 0; j < n; ++j)
  {
    for (std::size_t i = 0; i < n; ++i)
    {
      eigenvectors[i, j] = a_data[i + j * n];
    }
  }

  return {std::move(eigenvalues), std::move(eigenvectors), true};
#else
  return {vector<T>{}, matrix<T>{}, false};
#endif
}

} // namespace lam::linalg
