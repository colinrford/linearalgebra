
import std;
import lam.linearalgebra;

using namespace lam;

// ---- helpers ----------------------------------------------------------------

template<typename T>
T conj_val(T x)
{ return x; }
template<typename T>
std::complex<T> conj_val(std::complex<T> x)
{ return std::conj(x); }

void require(bool cond, const char* msg)
{
  if (!cond)
    throw std::runtime_error(msg);
}

// ---- individual checks -------------------------------------------------------

// Singular values, ordering, known values, and all svd_diagnostics fields
template<typename T>
void test_known_values()
{
  using R = decltype(std::abs(T{}));

  // A = [[1,2],[3,4],[5,6]]
  // A^T*A = [[35,44],[44,56]], trace=91, det=24
  // σ₁ = √((91+√(91²−4·24))/2) ≈ 9.52552,  σ₂ ≈ 0.51430
  // Also: σ₁²+σ₂² = ||A||_F² = 91
  matrix<T> A(3, 2);
  A[0, 0] = T{1};
  A[0, 1] = T{2};
  A[1, 0] = T{3};
  A[1, 1] = T{4};
  A[2, 0] = T{5};
  A[2, 1] = T{6};

  auto opt = linalg::svd(A);
  require(opt.has_value(), "svd returned nullopt");
  auto& res = *opt;

  require(res.s.size() == 2, "wrong singular value count");

  // Non-negative and descending
  require(res.s[0] >= res.s[1], "singular values not in descending order");
  require(res.s[1] >= R{0}, "negative singular value");

  // Known values — σ₁ = sqrt((91+sqrt(8185))/2), σ₂ = sqrt((91-sqrt(8185))/2)
  R kv_tol = std::is_same_v<R, float> ? R{1e-3} : R{1e-8};
  require(std::abs(res.s[0] - R{9.52551809}) < kv_tol, "s[0] wrong value");
  require(std::abs(res.s[1] - R{0.51430058}) < kv_tol, "s[1] wrong value");

  // Weyl check: σ₁²+σ₂² = ||A||_F²
  R sum_sq_tol = std::is_same_v<R, float> ? R{1e-3} : R{1e-9};
  require(std::abs(res.s[0] * res.s[0] + res.s[1] * res.s[1] - R{91}) < sum_sq_tol,
          "sum of squared singular values != ||A||_F^2");

  // Full diagnostics
  auto diag = linalg::svd_check(A, res);

  R recon_tol = std::is_same_v<R, float> ? R{1e-5} : R{1e-12};
  R orth_tol = std::is_same_v<R, float> ? R{1e-5} : R{1e-12};

  require(diag.reconstruction_error < recon_tol, "reconstruction error too large");
  require(diag.u_orthogonality_error < orth_tol, "U not orthogonal");
  require(diag.vt_orthogonality_error < orth_tol, "Vt not orthogonal");

  // Weyl bound should be small and positive
  require(diag.singular_value_bound > R{0}, "Weyl bound is zero");
  require(diag.condition_number >= R{1}, "condition number < 1");
  require(diag.min_singular_gap > R{0}, "min gap not positive for full-rank matrix");
}

// All four Moore-Penrose conditions
template<typename T>
void test_moore_penrose()
{
  using R = decltype(std::abs(T{}));

  matrix<T> A(3, 2);
  A[0, 0] = T{1};
  A[0, 1] = T{2};
  A[1, 0] = T{3};
  A[1, 1] = T{4};
  A[2, 0] = T{5};
  A[2, 1] = T{6};

  auto P = linalg::pinv(A); // 2×3

  R tol = std::is_same_v<R, float> ? R{1e-4} : R{1e-10};
  constexpr std::size_t m = 3, n = 2;

  // 1: A·pinv(A)·A = A
  {
    auto APA = A * P * A;
    for (std::size_t i = 0; i < m; ++i)
      for (std::size_t j = 0; j < n; ++j)
        require(std::abs(APA[i, j] - A[i, j]) < tol, "Moore-Penrose #1 failed");
  }

  // 2: pinv(A)·A·pinv(A) = pinv(A)
  {
    auto PAP = P * A * P;
    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = 0; j < m; ++j)
        require(std::abs(PAP[i, j] - P[i, j]) < tol, "Moore-Penrose #2 failed");
  }

  // 3: (A·pinv(A))^H = A·pinv(A)   — m×m Hermitian projection
  {
    auto AP = A * P;
    for (std::size_t i = 0; i < m; ++i)
      for (std::size_t j = 0; j < m; ++j)
        require(std::abs(conj_val(AP[j, i]) - AP[i, j]) < tol, "Moore-Penrose #3 failed");
  }

  // 4: (pinv(A)·A)^H = pinv(A)·A   — n×n Hermitian projection
  {
    auto PA = P * A;
    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = 0; j < n; ++j)
        require(std::abs(conj_val(PA[j, i]) - PA[i, j]) < tol, "Moore-Penrose #4 failed");
  }
}

// Rank-deficient matrix: threshold in matrix_rank and pinv should kick in
template<typename T>
void test_rank_deficient()
{
  using R = decltype(std::abs(T{}));

  // Rank-1 matrix: each row is a multiple of [1, 2]
  // A^T*A = 14·[[1,2],[2,4]], det=0, σ₁=√70≈8.366, σ₂=0 exactly
  matrix<T> A(3, 2);
  A[0, 0] = T{1};
  A[0, 1] = T{2};
  A[1, 0] = T{2};
  A[1, 1] = T{4};
  A[2, 0] = T{3};
  A[2, 1] = T{6};

  require(linalg::matrix_rank(A) == 1, "rank-1 matrix reported wrong rank");

  auto opt = linalg::svd(A);
  require(opt.has_value(), "svd of rank-deficient matrix failed");
  auto& res = *opt;

  // σ₂ should be numerically zero relative to σ₁
  R zero_tol = R{100} * std::numeric_limits<R>::epsilon() * res.s[0];
  require(res.s[1] < zero_tol, "s[1] not near zero for rank-1 matrix");

  // Diagnostics — reconstruction and orthogonality still hold
  auto diag = linalg::svd_check(A, res);
  R recon_tol = std::is_same_v<R, float> ? R{1e-5} : R{1e-12};
  R orth_tol = std::is_same_v<R, float> ? R{1e-5} : R{1e-12};
  require(diag.reconstruction_error < recon_tol, "reconstruction error too large (rank-deficient)");
  require(diag.u_orthogonality_error < orth_tol, "U not orthogonal (rank-deficient)");
  require(diag.vt_orthogonality_error < orth_tol, "Vt not orthogonal (rank-deficient)");
}

// lstsq: consistent system and normal-equations check for inconsistent system
template<typename T>
void test_lstsq()
{
  using R = decltype(std::abs(T{}));
  R tol = std::is_same_v<R, float> ? R{1e-4} : R{1e-10};

  // --- Consistent system: A*[1,2] = b exactly ---
  // A = [[1,0],[0,1],[1,1]], b = [1,2,3]
  {
    matrix<T> A(3, 2);
    A[0, 0] = T{1};
    A[0, 1] = T{0};
    A[1, 0] = T{0};
    A[1, 1] = T{1};
    A[2, 0] = T{1};
    A[2, 1] = T{1};

    vector<T> b(3);
    b[0] = T{1};
    b[1] = T{2};
    b[2] = T{3};

    auto x = linalg::lstsq(A, b);
    require(x.size() == 2, "lstsq output size wrong");
    require(std::abs(x[0] - T{1}) < tol, "lstsq x[0] wrong (consistent)");
    require(std::abs(x[1] - T{2}) < tol, "lstsq x[1] wrong (consistent)");

    // Residual should be near zero
    R res2 = R{0};
    for (std::size_t i = 0; i < 3; ++i)
    {
      T r = T{0};
      for (std::size_t j = 0; j < 2; ++j)
        r += A[i, j] * x[j];
      res2 += std::norm(r - b[i]);
    }
    require(std::sqrt(res2) < tol, "lstsq residual too large (consistent)");
  }

  // --- Inconsistent system: check normal equations A^T·(A·x - b) = 0 ---
  // A = [[1,0],[0,1],[1,1]], b = [2,3,4] — no exact solution
  {
    matrix<T> A(3, 2);
    A[0, 0] = T{1};
    A[0, 1] = T{0};
    A[1, 0] = T{0};
    A[1, 1] = T{1};
    A[2, 0] = T{1};
    A[2, 1] = T{1};

    vector<T> b(3);
    b[0] = T{2};
    b[1] = T{3};
    b[2] = T{4};

    auto x = linalg::lstsq(A, b);

    // Normal equations: A^T·A·x = A^T·b  <=>  sum_k conj(A[k,i])*(A[k,:]*x - b[k]) = 0 for each i
    R ne_tol = std::is_same_v<R, float> ? R{1e-4} : R{1e-10};
    for (std::size_t i = 0; i < 2; ++i)
    {
      T lhs = T{0};
      for (std::size_t k = 0; k < 3; ++k)
      {
        T Ax_k = T{0};
        for (std::size_t j = 0; j < 2; ++j)
          Ax_k += A[k, j] * x[j];
        lhs += conj_val(A[k, i]) * (Ax_k - b[k]);
      }
      require(std::abs(lhs) < ne_tol, "normal equations not satisfied (inconsistent)");
    }
  }
}

// ---- top-level ---------------------------------------------------------------

template<typename T>
void test_all()
{
  test_known_values<T>();
  test_moore_penrose<T>();
  test_rank_deficient<T>();
  test_lstsq<T>();
}

int main()
{
#ifdef LAM_USE_BLAS
  // The four type instantiations are fully independent; run them concurrently.
  auto f_f = std::async(std::launch::async, test_all<float>);
  auto f_d = std::async(std::launch::async, test_all<double>);
  auto f_cf = std::async(std::launch::async, test_all<std::complex<float>>);
  auto f_cd = std::async(std::launch::async, test_all<std::complex<double>>);
  f_f.get();
  f_d.get();
  f_cf.get();
  f_cd.get();
  std::println("SVD tests passed!");
#else
  std::println("BLAS not available, SVD tests skipped.");
#endif
  return 0;
}
