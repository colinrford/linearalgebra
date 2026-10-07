/*
 *  bm_eigen.cpp - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  lam vs Eigen on identical inputs. Each pair allocates its result the same way.
 *  lam's numbers depend on the build's BLAS setting; Eigen uses its own kernels.
 */

#include <Eigen/Dense>
#include <benchmark/benchmark.h>

import lam.linearalgebra;
import linalg.bench.support;

namespace lb = linalg::bench;

namespace
{
lam::linalg::vector<double> lam_vector(std::size_t n, std::uint64_t seed)
{
  lam::linalg::vector<double> v(n);
  lb::fill(v, n, seed);
  return v;
}

Eigen::VectorXd eigen_copy(const lam::linalg::vector<double>& v)
{
  Eigen::VectorXd e(static_cast<Eigen::Index>(v.size()));
  for (std::size_t i = 0; i < v.size(); ++i)
    e[static_cast<Eigen::Index>(i)] = v[i];
  return e;
}

template<lam::linalg::storage_layout Layout = lam::linalg::storage_layout::row_major>
lam::linalg::matrix<double, Layout> lam_matrix(std::size_t n, std::uint64_t seed)
{
  lam::linalg::matrix<double, Layout> m(n, n);
  lb::fill(m, n, n, seed);
  return m;
}

Eigen::MatrixXd eigen_copy(const lam::linalg::matrix<double>& m)
{
  Eigen::MatrixXd e(static_cast<Eigen::Index>(m.rows()), static_cast<Eigen::Index>(m.cols()));
  for (std::size_t i = 0; i < m.rows(); ++i)
    for (std::size_t j = 0; j < m.cols(); ++j)
      e(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) = m[i, j];
  return e;
}

void items(benchmark::State& state, std::size_t per_iteration)
{ state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(per_iteration)); }
} // namespace

static void bm_dot_lam(benchmark::State& state)
{
  const auto n = static_cast<std::size_t>(state.range(0));
  const auto a = lam_vector(n, 0xA);
  const auto b = lam_vector(n, 0xB);
  for (auto _ : state)
    benchmark::DoNotOptimize(lam::linalg::dot(a, b));
  items(state, n);
}

static void bm_dot_eigen(benchmark::State& state)
{
  const auto n = static_cast<std::size_t>(state.range(0));
  const auto a = eigen_copy(lam_vector(n, 0xA));
  const auto b = eigen_copy(lam_vector(n, 0xB));
  for (auto _ : state)
    benchmark::DoNotOptimize(a.dot(b));
  items(state, n);
}

static void bm_add_lam(benchmark::State& state)
{
  const auto n = static_cast<std::size_t>(state.range(0));
  const auto a = lam_vector(n, 0xA);
  const auto b = lam_vector(n, 0xB);
  for (auto _ : state)
  {
    auto c = a + b;
    benchmark::DoNotOptimize(c[0]);
    benchmark::ClobberMemory();
  }
  items(state, n);
}

static void bm_add_eigen(benchmark::State& state)
{
  const auto n = static_cast<std::size_t>(state.range(0));
  const auto a = eigen_copy(lam_vector(n, 0xA));
  const auto b = eigen_copy(lam_vector(n, 0xB));
  for (auto _ : state)
  {
    Eigen::VectorXd c = a + b;
    benchmark::DoNotOptimize(c[0]);
    benchmark::ClobberMemory();
  }
  items(state, n);
}

// Three operands: lam makes a temporary per operator, Eigen fuses into one pass.
static void bm_add3_lam(benchmark::State& state)
{
  const auto n = static_cast<std::size_t>(state.range(0));
  const auto a = lam_vector(n, 0xA);
  const auto b = lam_vector(n, 0xB);
  const auto c = lam_vector(n, 0xC);
  for (auto _ : state)
  {
    auto d = a + b + c;
    benchmark::DoNotOptimize(d[0]);
    benchmark::ClobberMemory();
  }
  items(state, n);
}

static void bm_add3_eigen(benchmark::State& state)
{
  const auto n = static_cast<std::size_t>(state.range(0));
  const auto a = eigen_copy(lam_vector(n, 0xA));
  const auto b = eigen_copy(lam_vector(n, 0xB));
  const auto c = eigen_copy(lam_vector(n, 0xC));
  for (auto _ : state)
  {
    Eigen::VectorXd d = a + b + c;
    benchmark::DoNotOptimize(d[0]);
    benchmark::ClobberMemory();
  }
  items(state, n);
}

template<lam::linalg::storage_layout Layout>
static void bm_gemm_lam(benchmark::State& state)
{
  const auto n = static_cast<std::size_t>(state.range(0));
  const auto a = lam_matrix<Layout>(n, 0xA);
  const auto b = lam_matrix<Layout>(n, 0xB);
  for (auto _ : state)
  {
    auto c = a * b;
    benchmark::DoNotOptimize(c[0, 0]);
    benchmark::ClobberMemory();
  }
  items(state, 2 * n * n * n);
}

static void bm_gemm_eigen(benchmark::State& state)
{
  const auto n = static_cast<std::size_t>(state.range(0));
  const auto a = eigen_copy(lam_matrix(n, 0xA));
  const auto b = eigen_copy(lam_matrix(n, 0xB));
  for (auto _ : state)
  {
    Eigen::MatrixXd c = a * b;
    benchmark::DoNotOptimize(c(0, 0));
    benchmark::ClobberMemory();
  }
  items(state, 2 * n * n * n);
}

BENCHMARK(bm_dot_lam)->RangeMultiplier(8)->Range(64, 1 << 18);
BENCHMARK(bm_dot_eigen)->RangeMultiplier(8)->Range(64, 1 << 18);
BENCHMARK(bm_add_lam)->RangeMultiplier(8)->Range(64, 1 << 18);
BENCHMARK(bm_add_eigen)->RangeMultiplier(8)->Range(64, 1 << 18);
BENCHMARK(bm_add3_lam)->RangeMultiplier(8)->Range(64, 1 << 18);
BENCHMARK(bm_add3_eigen)->RangeMultiplier(8)->Range(64, 1 << 18);
BENCHMARK_TEMPLATE(bm_gemm_lam, lam::linalg::storage_layout::row_major)->RangeMultiplier(2)->Range(16, 512);
BENCHMARK_TEMPLATE(bm_gemm_lam, lam::linalg::storage_layout::col_major)->RangeMultiplier(2)->Range(16, 512);
BENCHMARK(bm_gemm_eigen)->RangeMultiplier(2)->Range(16, 512);
