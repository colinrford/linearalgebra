/*
 *  bm_allocator_profiles.cpp - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  The allocator-by-profile centerpiece: ONE lam::linalg::vector type, the
 *  allocator as the variable, across three application allocation profiles.
 *  The matrix/vector code never changes -- only which allocator is plugged in,
 *  and the right choice is dictated by the workload's shape and lifetime.
 *
 *    games      -- many small, transient objects; arena + per-frame reset
 *    graphics3d -- many small objects, allocate-then-free churn; pool recycling
 *    ml         -- few large, long-lived objects; std::allocator is already fine
 */

#include <benchmark/benchmark.h>

import std;
import lam.linearalgebra;
import linalg.bench.support;

namespace lb = linalg::bench;

using pmr_vec = lam::linalg::vector<double, std::pmr::polymorphic_allocator<double>>;
using heap_vec = lam::linalg::vector<double>;

// A flattened 4x4 transform's worth of doubles -- the small-object size that
// dominates games/3D allocation traffic.
inline constexpr std::size_t transform_dim = 16;

// ---------------------------------------------------------------------------
// games profile: `batch` small transient objects per frame, thrown away each
// frame. Baseline churns the general heap; the arena bumps a pointer and the
// per-frame reset() rewinds it in O(1). reset() is INSIDE the timed loop --
// it's part of the workload, and the arena still wins despite paying it.
// ---------------------------------------------------------------------------
static void bm_games_heap(benchmark::State& state)
{
  const std::size_t batch = static_cast<std::size_t>(state.range(0));
  for (auto _ : state)
    for (std::size_t k = 0; k < batch; ++k)
    {
      heap_vec v(transform_dim);
      lb::fill(v, transform_dim, k);
      benchmark::DoNotOptimize(v[0]);
    }
  state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(batch));
}

static void bm_games_arena(benchmark::State& state)
{
  const std::size_t batch = static_cast<std::size_t>(state.range(0));
  // Buffer sized to hold a full frame's batch (with slack) so a frame never
  // spills to the upstream heap; reset() reclaims it wholesale each frame.
  lb::monotonic_arena<double> arena{batch * transform_dim * sizeof(double) * 2 + 4096};

  for (auto _ : state)
  {
    for (std::size_t k = 0; k < batch; ++k)
    {
      pmr_vec v(transform_dim, arena.allocator());
      lb::fill(v, transform_dim, k);
      benchmark::DoNotOptimize(v[0]);
    }
    arena.reset(); // the per-frame rewind
  }
  state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(batch));
}

// ---------------------------------------------------------------------------
// graphics3d profile: same small-object traffic, but allocate-then-free churn
// (objects freed within the frame). A pool resource recycles same-sized blocks
// instead of round-tripping through malloc/free. Pool persists across frames.
// ---------------------------------------------------------------------------
static void bm_graphics_heap(benchmark::State& state)
{
  const std::size_t batch = static_cast<std::size_t>(state.range(0));
  for (auto _ : state)
    for (std::size_t k = 0; k < batch; ++k)
    {
      heap_vec v(transform_dim);
      lb::fill(v, transform_dim, k);
      benchmark::DoNotOptimize(v[0]);
    }
  state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(batch));
}

static void bm_graphics_pool(benchmark::State& state)
{
  const std::size_t batch = static_cast<std::size_t>(state.range(0));
  std::pmr::unsynchronized_pool_resource pool;
  std::pmr::polymorphic_allocator<double> alloc{&pool};

  for (auto _ : state)
    for (std::size_t k = 0; k < batch; ++k)
    {
      pmr_vec v(transform_dim, alloc); // freed at scope exit -> pool recycles
      lb::fill(v, transform_dim, k);
      benchmark::DoNotOptimize(v[0]);
    }
  state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(batch));
}

// ---------------------------------------------------------------------------
// ml profile: few large, long-lived operands; each op allocates one big result.
// Allocation is a vanishing fraction of the compute + copy, so the general heap
// is already fine -- the honest, useful "and here the allocator stops mattering"
// data point. Size axis is the vector length (contrast: games/3D swept batch).
// ---------------------------------------------------------------------------
static void bm_ml_heap(benchmark::State& state)
{
  const std::size_t dim = static_cast<std::size_t>(state.range(0));
  heap_vec a(dim), b(dim);
  lb::fill(a, dim, 1);
  lb::fill(b, dim, 2);

  for (auto _ : state)
  {
    auto c = a + b; // one large allocation per op
    benchmark::DoNotOptimize(c[0]);
    benchmark::ClobberMemory();
  }
  state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(dim));
}

// games / graphics3d: sweep objects-per-frame (small fixed object size).
BENCHMARK(bm_games_heap)->RangeMultiplier(8)->Range(8, 512);
BENCHMARK(bm_games_arena)->RangeMultiplier(8)->Range(8, 512);
BENCHMARK(bm_graphics_heap)->RangeMultiplier(8)->Range(8, 512);
BENCHMARK(bm_graphics_pool)->RangeMultiplier(8)->Range(8, 512);
// ml: sweep operand length (few, large, long-lived).
BENCHMARK(bm_ml_heap)->RangeMultiplier(8)->Range(1024, 262144);
