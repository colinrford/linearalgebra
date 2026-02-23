/*
 *  benchmark_svd.cpp - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  Measures gesvd (QR) vs gesdd (divide-and-conquer) across a range of sizes
 *  to find where the D&C crossover makes sense on this machine.
 *
 *  Results — 2018 MacBook Pro 13", Intel Core i7-8559U 2.7 GHz, Apple Accelerate:
 *    Crossover at min(m,n) = 30:
 *      n=25 → gesvd wins (ratio 0.97)
 *      n=30 → gesdd wins (ratio 1.24), advantage grows to 13x+ at n=2000
 *    Rectangular matrices show the same threshold (based on min dimension).
 *    Recommended threshold for *gesdd_ dispatch: min(m,n) >= 30
 */

#include <Accelerate/Accelerate.h>

import std;

// ---- LAPACK declarations ----------------------------------------------------

extern "C"
{
  // QR-based SVD
  void dgesvd_(const char* jobu, const char* jobvt, const int* m, const int* n, double* a, const int* lda, double* s,
               double* u, const int* ldu, double* vt, const int* ldvt, double* work, const int* lwork, int* info);

  // Divide-and-conquer SVD
  void dgesdd_(const char* jobz, const int* m, const int* n, double* a, const int* lda, double* s, double* u,
               const int* ldu, double* vt, const int* ldvt, double* work, const int* lwork, int* iwork, int* info);
}

// ---- timing helper ----------------------------------------------------------

// Returns the minimum wall time in microseconds over `reps` calls to f().
// Matrix data is re-copied from `src` before each call so each run sees a fresh input.
template<typename Setup, typename Func>
double min_us(int reps, Setup&& setup, Func&& f)
{
  double best = std::numeric_limits<double>::max();
  for (int rep = 0; rep < reps; ++rep)
  {
    setup();
    auto t0 = std::chrono::steady_clock::now();
    f();
    auto t1 = std::chrono::steady_clock::now();
    double us = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count()) / 1e3;
    best = std::min(best, us);
  }
  return best;
}

// ---- per-size benchmark -----------------------------------------------------

struct result
{
  int m, n;
  double gesvd_us;
  double gesdd_us;
};

result bench_size(int m, int n, int reps)
{
  const int k = std::min(m, n);
  const int mn = m * n;

  // Source matrix filled once with deterministic values
  std::vector<double> src(static_cast<std::size_t>(mn));
  {
    std::mt19937_64 rng{0xDEADBEEF};
    std::uniform_real_distribution<double> dist{-1.0, 1.0};
    for (auto& v : src)
      v = dist(rng);
  }

  // Working buffers (sized for the largest possible workspace query)
  std::vector<double> a(static_cast<std::size_t>(mn));
  std::vector<double> s(static_cast<std::size_t>(k));
  std::vector<double> u(static_cast<std::size_t>(m * m));
  std::vector<double> vt(static_cast<std::size_t>(n * n));

  // ---- gesvd ----------------------------------------------------------------
  double gesvd_us = 0;
  {
    // Workspace query
    char jobu = 'A', jobvt = 'A';
    int lda = m, ldu = m, ldvt = n, lwork = -1, info = 0;
    double work_query = 0;
    std::copy(src.begin(), src.end(), a.begin());
    dgesvd_(&jobu, &jobvt, &m, &n, a.data(), &lda, s.data(), u.data(), &ldu, vt.data(), &ldvt, &work_query, &lwork,
            &info);
    lwork = static_cast<int>(work_query);
    std::vector<double> work(static_cast<std::size_t>(lwork));

    gesvd_us = min_us(
      reps, [&]() { std::copy(src.begin(), src.end(), a.begin()); },
      [&]() {
        dgesvd_(&jobu, &jobvt, &m, &n, a.data(), &lda, s.data(), u.data(), &ldu, vt.data(), &ldvt, work.data(), &lwork,
                &info);
      });
  }

  // ---- gesdd ----------------------------------------------------------------
  double gesdd_us = 0;
  {
    // Workspace query
    char jobz = 'A';
    int lda = m, ldu = m, ldvt = n, lwork = -1, info = 0;
    double work_query = 0;
    std::vector<int> iwork(static_cast<std::size_t>(8 * k));
    std::copy(src.begin(), src.end(), a.begin());
    dgesdd_(&jobz, &m, &n, a.data(), &lda, s.data(), u.data(), &ldu, vt.data(), &ldvt, &work_query, &lwork,
            iwork.data(), &info);
    lwork = static_cast<int>(work_query);
    std::vector<double> work(static_cast<std::size_t>(lwork));

    gesdd_us = min_us(
      reps, [&]() { std::copy(src.begin(), src.end(), a.begin()); },
      [&]() {
        dgesdd_(&jobz, &m, &n, a.data(), &lda, s.data(), u.data(), &ldu, vt.data(), &ldvt, work.data(), &lwork,
                iwork.data(), &info);
      });
  }

  return {m, n, gesvd_us, gesdd_us};
}

// ---- driver -----------------------------------------------------------------

int main()
{
  constexpr int reps = 7;

  // Square sizes — dense near the expected crossover (25–75) for precise threshold
  std::vector<int> sq_sizes = {10, 20,  25,  30,  35,  40,  45,  50,  55,   60,   70,  75,
                               80, 100, 150, 200, 300, 400, 500, 750, 1000, 1500, 2000};

  // Rectangular sizes (tall 2:1 and wide 1:2)
  std::vector<std::pair<int, int>> rect_sizes = {
    {200, 100}, {400, 200}, {800, 400}, {1000, 500}, {100, 200}, {200, 400}, {400, 800}, {500, 1000},
  };

  std::println("gesvd vs gesdd — double, full singular vectors, {} reps (min time)\n", reps);
  std::println("{:>6}  {:>6}  {:>10}  {:>10}  {:>6}  {}", "m", "n", "gesvd µs", "gesdd µs", "ratio", "winner");
  std::println("{}", std::string(60, '-'));

  auto print_row = [](result r) {
    double ratio = r.gesvd_us / r.gesdd_us;
    const char* winner = (r.gesdd_us < r.gesvd_us) ? "gesdd ✓" : "gesvd";
    std::println("{:>6}  {:>6}  {:>10.1f}  {:>10.1f}  {:>6.2f}  {}", r.m, r.n, r.gesvd_us, r.gesdd_us, ratio, winner);
  };

  std::println("-- square --");
  for (int sz : sq_sizes)
    print_row(bench_size(sz, sz, reps));

  std::println("\n-- tall (2:1) --");
  for (auto [m, n] : rect_sizes)
    if (m > n)
      print_row(bench_size(m, n, reps));

  std::println("\n-- wide (1:2) --");
  for (auto [m, n] : rect_sizes)
    if (m < n)
      print_row(bench_size(m, n, reps));

  return 0;
}
