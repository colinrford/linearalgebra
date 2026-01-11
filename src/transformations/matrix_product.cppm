export module lam.linearalgebra:matrix.product;

import std;
import :config;
import :matrix;

namespace lam::linalg
{

#ifdef LAM_USE_BLAS
// Manually declare CBLAS types/constants at namespace scope
enum CBLAS_ORDER
{
  CblasRowMajor = 101,
  CblasColMajor = 102
};
enum CBLAS_TRANSPOSE
{
  CblasNoTrans = 111,
  CblasTrans = 112,
  CblasConjTrans = 113
};

extern "C"
{
  void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB,
                   const int M, const int N, const int K, const double alpha, const double* A, const int lda,
                   const double* B, const int ldb, const double beta, double* C, const int ldc);
  void cblas_sgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB,
                   const int M, const int N, const int K, const float alpha, const float* A, const int lda,
                   const float* B, const int ldb, const float beta, float* C, const int ldc);
}

// Type-safe dispatch helper
template<typename T>
struct blas_dispatcher;

template<>
struct blas_dispatcher<double>
{
  static void gemm(CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA, CBLAS_TRANSPOSE TransB, int M, int N, int K, double alpha,
                   const double* A, int lda, const double* B, int ldb, double beta, double* C, int ldc)
  {
    cblas_dgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
  }
};

template<>
struct blas_dispatcher<float>
{
  static void gemm(CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA, CBLAS_TRANSPOSE TransB, int M, int N, int K, float alpha,
                   const float* A, int lda, const float* B, int ldb, float beta, float* C, int ldc)
  {
    cblas_sgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
  }
};
#endif


// Helper for generic multiplication (constexpr-compatible)
template<typename T, typename Alloc, storage_layout LayoutA, storage_layout LayoutB>
constexpr void generic_matrix_multiply(const matrix<T, Alloc, LayoutA>& a, const matrix<T, Alloc, LayoutB>& b,
                                       matrix<T, Alloc, storage_layout::row_major>& res)
{
  // Case 1: RowMajor * ColMajor -> Inner Product Loop (i, j, k)
  if constexpr (LayoutA == storage_layout::row_major && LayoutB == storage_layout::col_major)
  {
    for (std::size_t i = 0; i < a.rows(); ++i)
    {
      for (std::size_t j = 0; j < b.cols(); ++j)
      {
        T sum = T{0};
        for (std::size_t k = 0; k < a.cols(); ++k)
        {
          sum += a[i, k] * b[k, j];
        }
        res[i, j] = sum;
      }
    }
  }
  // Case 2: ColMajor * ColMajor -> Inner Product Loop (i, j, k)
  else if constexpr (LayoutA == storage_layout::col_major && LayoutB == storage_layout::col_major)
  {
    for (std::size_t i = 0; i < a.rows(); ++i)
    {
      for (std::size_t j = 0; j < b.cols(); ++j)
      {
        T sum = T{0};
        for (std::size_t k = 0; k < a.cols(); ++k)
        {
          sum += a[i, k] * b[k, j];
        }
        res[i, j] = sum;
      }
    }
  }
  // Case 3: Cache-Tiled Multiplication
  else
  {
    constexpr std::size_t BLOCK = 64;
    const std::size_t M = a.rows();
    const std::size_t N = b.cols();
    const std::size_t K = a.cols();

    for (std::size_t i = 0; i < M; ++i)
      for (std::size_t j = 0; j < N; ++j)
        res[i, j] = T{0};

    // Tiled loops
    for (std::size_t ii = 0; ii < M; ii += BLOCK)
    {
      const std::size_t i_end = std::min(ii + BLOCK, M);
      for (std::size_t kk = 0; kk < K; kk += BLOCK)
      {
        const std::size_t k_end = std::min(kk + BLOCK, K);
        for (std::size_t jj = 0; jj < N; jj += BLOCK)
        {
          const std::size_t j_end = std::min(jj + BLOCK, N);

          for (std::size_t i = ii; i < i_end; ++i)
          {
            for (std::size_t k = kk; k < k_end; ++k)
            {
              T r = a[i, k];
              for (std::size_t j = jj; j < j_end; ++j)
              {
                res[i, j] += r * b[k, j];
              }
            }
          }
        }
      }
    }
  }
}

export template<typename T, typename Alloc, storage_layout LayoutA, storage_layout LayoutB>
constexpr matrix<T, Alloc, storage_layout::row_major> operator*(const matrix<T, Alloc, LayoutA>& a,
                                                                const matrix<T, Alloc, LayoutB>& b)
{
  if (a.cols() != b.rows())
    throw matrix_exception::dim_mismatch();

  matrix<T, Alloc, storage_layout::row_major> res(a.rows(), b.cols());

  if consteval
  {
    generic_matrix_multiply(a, b, res);
  }
  else
  {
    // Hybrid BLAS Config
    if constexpr (config::use_blas && (std::is_same_v<T, double> || std::is_same_v<T, float>))
    {
#ifdef LAM_USE_BLAS
      const auto transA = (LayoutA == storage_layout::row_major) ? CblasNoTrans : CblasTrans;
      const auto transB = (LayoutB == storage_layout::row_major) ? CblasNoTrans : CblasTrans;
      // Dimensions
      const int M = static_cast<int>(a.rows());
      const int N = static_cast<int>(b.cols());
      const int K = static_cast<int>(a.cols());
      // Strides
      const int lda = (LayoutA == storage_layout::row_major) ? K : M;
      const int ldb = (LayoutB == storage_layout::row_major) ? N : K;
      const int ldc = N;

      blas_dispatcher<T>::gemm(CblasRowMajor, transA, transB, M, N, K, T{1.0}, a.begin(), lda, b.begin(), ldb, T{0.0},
                               res.begin(), ldc);
      return res;
#else
      generic_matrix_multiply(a, b, res);
#endif
    }
    else
    {
      generic_matrix_multiply(a, b, res);
    }
  }

  return res;
}

} // namespace lam::linalg
