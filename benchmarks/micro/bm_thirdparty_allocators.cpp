/*
 *  bm_thirdparty_allocators.cpp - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  The allocator-by-profile benchmark's companion: same workload, but the
 *  competitors are other people's allocators rather than std's. One
 *  lam::linalg::vector type throughout; only the allocator changes.
 *
 *  The workload is the graphics3d churn profile from bm_allocator_profiles --
 *  many small objects allocated and freed within a frame -- because that is
 *  where a general-purpose allocator has something to prove against a pool.
 *  The arena is not here: it wins by not freeing, which is a different claim.
 *
 *  Both dependencies are optional and quarantined; the parent CMakeLists only
 *  builds this when find_package found them.
 */

#include <benchmark/benchmark.h>
#include <boost/container/pmr/memory_resource.hpp>
#include <boost/container/pmr/unsynchronized_pool_resource.hpp>
#include <boost/pool/pool_alloc.hpp>
#include <oneapi/tbb/cache_aligned_allocator.h>
#include <oneapi/tbb/scalable_allocator.h>

import std;
import lam.linearalgebra;
import linalg.bench.support;

namespace lb = linalg::bench;

// Same small-object size as bm_allocator_profiles: a flattened 4x4 transform.
inline constexpr std::size_t transform_dim = 16;

// Boost's memory_resource is a separate base from std's, so reaching it takes a local adaptor.
class boost_resource_adaptor final : public std::pmr::memory_resource
{
public:
  explicit boost_resource_adaptor(boost::container::pmr::memory_resource& upstream) : m_upstream{&upstream} {}

private:
  void* do_allocate(std::size_t bytes, std::size_t alignment) override
  { return m_upstream->allocate(bytes, alignment); }
  void do_deallocate(void* p, std::size_t bytes, std::size_t alignment) override
  { m_upstream->deallocate(p, bytes, alignment); }
  bool do_is_equal(const std::pmr::memory_resource& other) const noexcept override { return this == &other; }

  boost::container::pmr::memory_resource* m_upstream;
};

// The churn loop, identical for every allocator: allocate `batch` small vectors, fill, destroy.
template<typename Vec, typename MakeAlloc>
static void churn(benchmark::State& state, MakeAlloc make_alloc)
{
  const std::size_t batch = static_cast<std::size_t>(state.range(0));
  auto alloc = make_alloc();

  for (auto _ : state)
    for (std::size_t k = 0; k < batch; ++k)
    {
      Vec v(transform_dim, alloc);
      lb::fill(v, transform_dim, k);
      benchmark::DoNotOptimize(v[0]);
    }
  state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(batch));
}

using pmr_vec = lam::linalg::vector<double, std::pmr::polymorphic_allocator<double>>;
using tbb_scalable_vec = lam::linalg::vector<double, tbb::scalable_allocator<double>>;
using tbb_aligned_vec = lam::linalg::vector<double, tbb::cache_aligned_allocator<double>>;
using boost_pool_vec = lam::linalg::vector<double, boost::pool_allocator<double>>;

// Baseline: the general heap, via std::allocator. Every row below is measured against this.
static void bm_churn_heap(benchmark::State& state)
{
  const std::size_t batch = static_cast<std::size_t>(state.range(0));
  for (auto _ : state)
    for (std::size_t k = 0; k < batch; ++k)
    {
      lam::linalg::vector<double> v(transform_dim);
      lb::fill(v, transform_dim, k);
      benchmark::DoNotOptimize(v[0]);
    }
  state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(batch));
}

// std's own pool, the in-tree answer to this profile.
static void bm_churn_std_pool(benchmark::State& state)
{
  std::pmr::unsynchronized_pool_resource pool;
  churn<pmr_vec>(state, [&] { return std::pmr::polymorphic_allocator<double>{&pool}; });
}

// Boost's pool through the adaptor; the adaptor's cost is included on purpose.
static void bm_churn_boost_pool(benchmark::State& state)
{
  boost::container::pmr::unsynchronized_pool_resource pool;
  boost_resource_adaptor adaptor(pool);
  churn<pmr_vec>(state, [&] { return std::pmr::polymorphic_allocator<double>{&adaptor}; });
}

// Boost's plain Allocator: no adaptor, no type erasure, singleton pool behind it.
static void bm_churn_boost_pool_allocator(benchmark::State& state)
{
  churn<boost_pool_vec>(state, [] { return boost::pool_allocator<double>{}; });
}

// TBB's scalable allocator: a malloc replacement, not specialized to this size class.
static void bm_churn_tbb_scalable(benchmark::State& state)
{
  churn<tbb_scalable_vec>(state, [] { return tbb::scalable_allocator<double>{}; });
}

// Pays for 64-byte alignment on every block; here to show what that costs at this size.
static void bm_churn_tbb_cache_aligned(benchmark::State& state)
{
  churn<tbb_aligned_vec>(state, [] { return tbb::cache_aligned_allocator<double>{}; });
}

// Same sweep as the graphics3d rows in bm_allocator_profiles, so the two are comparable.
BENCHMARK(bm_churn_heap)->RangeMultiplier(8)->Range(8, 512);
BENCHMARK(bm_churn_std_pool)->RangeMultiplier(8)->Range(8, 512);
BENCHMARK(bm_churn_boost_pool)->RangeMultiplier(8)->Range(8, 512);
BENCHMARK(bm_churn_boost_pool_allocator)->RangeMultiplier(8)->Range(8, 512);
BENCHMARK(bm_churn_tbb_scalable)->RangeMultiplier(8)->Range(8, 512);
BENCHMARK(bm_churn_tbb_cache_aligned)->RangeMultiplier(8)->Range(8, 512);
