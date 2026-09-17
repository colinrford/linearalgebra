/*
 *  bm_smoke.cpp - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  Smoke benchmark: proves the micro/ plumbing works end to end --
 *  import lam.linearalgebra + import linalg.bench.support + link google/benchmark.
 *  The DoNotOptimize barrier lives here (in the gb-including TU), not in the module.
 */

#include <benchmark/benchmark.h>

import lam.linearalgebra;
import linalg.bench.support;

namespace lb = linalg::bench;

// Heap-allocated vector add, deterministic fill, over a size sweep. Just enough
// to exercise module import, the support module, and the gb link/counters.
static void bm_vector_add_heap(benchmark::State& state)
{
  const std::size_t n = static_cast<std::size_t>(state.range(0));

  lam::linalg::vector<double> a(n);
  lam::linalg::vector<double> b(n);
  lb::fill(a, n, 0xA);
  lb::fill(b, n, 0xB);

  for (auto _ : state)
  {
    auto c = a + b;
    benchmark::DoNotOptimize(c[0]);
    benchmark::ClobberMemory();
  }
  state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(n));
}
BENCHMARK(bm_vector_add_heap)->RangeMultiplier(8)->Range(8, 4096);
