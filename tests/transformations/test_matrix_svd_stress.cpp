
import std;
import lam.linearalgebra;

using namespace lam;

struct rng_state
{
  std::mt19937 engine{0xDEADBEEF};
  std::uniform_real_distribution<double> dist{-10.0, 10.0};
};

template<typename T>
struct value_gen;

template<>
struct value_gen<float>
{
  static float get(rng_state& rng) { return static_cast<float>(rng.dist(rng.engine)); }
};

template<>
struct value_gen<double>
{
  static double get(rng_state& rng) { return rng.dist(rng.engine); }
};

template<>
struct value_gen<std::complex<float>>
{
  static std::complex<float> get(rng_state& rng)
  { return {static_cast<float>(rng.dist(rng.engine)), static_cast<float>(rng.dist(rng.engine))}; }
};

template<>
struct value_gen<std::complex<double>>
{
  static std::complex<double> get(rng_state& rng) { return {rng.dist(rng.engine), rng.dist(rng.engine)}; }
};

// ---- per-case checks --------------------------------------------------------

template<typename T>
void run_stress_test(std::size_t rows, std::size_t cols)
{
  using R = decltype(std::abs(T{}));

  // Deterministic seed derived from dimensions — independent of task ordering
  rng_state rng;
  rng.engine.seed(0xDEADBEEF ^ (static_cast<std::uint32_t>(rows) * 2654435761u + static_cast<std::uint32_t>(cols)));

  matrix<T> A(rows, cols);
  for (std::size_t i = 0; i < rows; ++i)
    for (std::size_t j = 0; j < cols; ++j)
      A[i, j] = value_gen<T>::get(rng);

  auto opt = linalg::svd(A);
  if (!opt)
    throw std::runtime_error("SVD failed entirely during stress test.");
  auto& res = *opt;

  const std::size_t k = res.s.size();

  // Singular values: non-negative and descending
  for (std::size_t i = 0; i < k; ++i)
    if (res.s[i] < R{0})
      throw std::runtime_error("Negative singular value in stress test.");
  for (std::size_t i = 0; i + 1 < k; ++i)
    if (res.s[i] < res.s[i + 1] - R{100} * std::numeric_limits<R>::epsilon() * res.s[i])
      throw std::runtime_error("Singular values not descending in stress test.");

  // Reconstruction error (O(mn·k) — always fast)
  {
    R frob_A_sq = R{0}, frob_res_sq = R{0};
    for (std::size_t i = 0; i < rows; ++i)
      for (std::size_t j = 0; j < cols; ++j)
      {
        frob_A_sq += std::norm(A[i, j]);
        T rec = T{0};
        for (std::size_t r = 0; r < k; ++r)
          rec += res.u[i, r] * static_cast<T>(res.s[r]) * res.vt[r, j];
        frob_res_sq += std::norm(A[i, j] - rec);
      }
    R rel_err = (frob_A_sq > R{0}) ? std::sqrt(frob_res_sq / frob_A_sq) : std::sqrt(frob_res_sq);

    // Expected backward error: O(min(rows,cols)·ε)
    R recon_tol = std::is_same_v<R, float> ? static_cast<R>(std::min(rows, cols)) * R{1e-4}
                                           : static_cast<R>(std::min(rows, cols)) * R{1e-13};

    if (rel_err > recon_tol)
    {
      std::println("  Reconstruction error {} > {} for {}x{}", rel_err, recon_tol, rows, cols);
      throw std::runtime_error("Stress test: reconstruction error too large");
    }
  }

  // Orthogonality checks via svd_check (O(max(m,n)³) — only for smaller matrices)
  if (std::max(rows, cols) <= 100)
  {
    auto diag = linalg::svd_check(A, res);
    R orth_tol = std::is_same_v<R, float> ? R{std::sqrt(static_cast<R>(std::max(rows, cols)))} * R{1e-4}
                                          : R{std::sqrt(static_cast<R>(std::max(rows, cols)))} * R{1e-12};

    if (diag.u_orthogonality_error > orth_tol)
    {
      std::println("  U orthogonality error {} > {} for {}x{}", diag.u_orthogonality_error, orth_tol, rows, cols);
      throw std::runtime_error("Stress test: U orthogonality error too large");
    }
    if (diag.vt_orthogonality_error > orth_tol)
    {
      std::println("  Vt orthogonality error {} > {} for {}x{}", diag.vt_orthogonality_error, orth_tol, rows, cols);
      throw std::runtime_error("Stress test: Vt orthogonality error too large");
    }
  }

  // Pseudoinverse property: A·pinv(A)·A = A  (existing check, retained)
  {
    auto pinv_A = linalg::pinv(A);
    auto A_PA_A = A * pinv_A * A;

    R max_err = R{0};
    for (std::size_t i = 0; i < rows; ++i)
      for (std::size_t j = 0; j < cols; ++j)
        max_err = std::max(max_err, std::abs(A_PA_A[i, j] - A[i, j]));

    R pinv_tol = std::is_same_v<R, float> ? R{2e-2} : R{1e-9};
    if (max_err > pinv_tol)
    {
      std::println("  Pinv error {} > {} for {}x{}", max_err, pinv_tol, rows, cols);
      throw std::runtime_error("Stress test: pseudoinverse check failed");
    }
  }
}

// ---- driver -----------------------------------------------------------------

template<typename T>
void stress_test_type(std::string_view type_name)
{
  // Square and rectangular (tall and wide); sizes up to 500
  std::vector<std::pair<std::size_t, std::size_t>> sizes = {{10, 10},   {50, 50},  {100, 100}, {200, 200},
                                                            {500, 500}, {20, 10},  {100, 50},  {500, 100},
                                                            {10, 20},   {50, 100}, {100, 500}};

  std::vector<std::future<void>> futs;
  futs.reserve(sizes.size());

  for (auto [rows, cols] : sizes)
  {
    // Skip the very largest cases for complex types (16-byte scalars)
    if (rows * cols > 100000 && sizeof(T) > 8)
      continue;
    std::println("  [{}] {}x{}...", type_name, rows, cols);
    futs.push_back(std::async(std::launch::async, [rows, cols]() { run_stress_test<T>(rows, cols); }));
  }

  for (auto& f : futs)
    f.get();
}

int main()
{
#ifdef LAM_USE_BLAS
  try
  {
    // Run all four type families in parallel
    auto f_f = std::async(std::launch::async, stress_test_type<float>, "float");
    auto f_d = std::async(std::launch::async, stress_test_type<double>, "double");
    auto f_cf = std::async(std::launch::async, stress_test_type<std::complex<float>>, "complex<float>");
    auto f_cd = std::async(std::launch::async, stress_test_type<std::complex<double>>, "complex<double>");
    f_f.get();
    f_d.get();
    f_cf.get();
    f_cd.get();
    std::println("SVD stress tests passed!");
  }
  catch (const std::exception& e)
  {
    std::println(std::cerr, "Stress test failed: {}", e.what());
    return 1;
  }
#else
  std::println("BLAS not available, SVD stress tests skipped.");
#endif
  return 0;
}
