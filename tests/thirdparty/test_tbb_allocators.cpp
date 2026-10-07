/*
 *  test_tbb_allocators.cpp - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  test_tbb_allocators is a c++ module
 *
 *  Correctness only: does lam.linalg behave when handed an allocator it has
 *  never seen? Nothing here measures speed. TBB is optional -- this whole
 *  directory is skipped when find_package(TBB) fails.
 */

#include <oneapi/tbb/cache_aligned_allocator.h>
#include <oneapi/tbb/scalable_allocator.h>

import std;
import lam.linearalgebra;

import linalg.test.pmr_suite;

using lam::linalg::matrix;
using lam::linalg::storage_layout;

int failures = 0;

void check(bool ok, std::string_view what)
{
  if (!ok)
  {
    std::println("FAILED: {}", what);
    ++failures;
  }
}

// What lam.linalg assumes of any allocator; a third-party change should fail here, not at runtime.
template<typename A>
void allocator_shape(std::string_view name)
{
  using traits = std::allocator_traits<A>;
  static_assert(std::same_as<typename traits::value_type, double>);
  static_assert(std::same_as<typename traits::pointer, double*>, "lam.linalg stores raw T*, not a fancy pointer");
  static_assert(std::same_as<typename traits::template rebind_alloc<int>::value_type, int>);
  static_assert(std::copy_constructible<A>);
  static_assert(std::equality_comparable<A>);
  std::println("  {}: always_equal={}, pocma={}", name, traits::is_always_equal::value,
               traits::propagate_on_container_move_assignment::value);
}

template<typename A>
void container_suite(std::string_view name)
{
  using alloc = A;
  using vec = lam::linalg::vector<double, alloc>;
  using mat = matrix<double, storage_layout::row_major, alloc>;

  vec v({3.0, 4.0, 0.0}, alloc{});
  check(v.size() == 3, "vector size");
  check(v.norm() == 5.0, "vector norm");

  vec w = v + v;
  check(w[1] == 8.0, "vector operator+");
  check(w.get_allocator() == v.get_allocator(), "result kept the operand's allocator");

  // Move-assign across two instances of a stateless allocator takes the steal path.
  vec z(alloc{});
  z = std::move(w);
  check(z[1] == 8.0, "vector move-assign");
  check(w.size() == 0, "moved-from vector is empty");

  mat m(3, 3, alloc{});
  for (std::size_t i = 0; i < 3; ++i)
    m[i, i] = 2.0;
  auto p = m * m;
  check(p[1, 1] == 4.0, "matrix product");
  check(p.get_allocator() == m.get_allocator(), "matrix result kept the operand's allocator");

  auto mv = m * v;
  check(mv[0] == 6.0, "matrix * vector");

  // Decompositions rebind this allocator for their own workspaces.
  auto s = lam::linalg::svd(m);
  check(s.has_value(), "svd succeeded");
  if (s)
    check(std::abs(s->s[0] - 2.0) < 1e-12, "svd singular value");

  auto e = lam::linalg::symmetric_eigen(m);
  check(e.success, "symmetric_eigen succeeded");
  if (e.success)
    check(std::abs(e.eigenvalues[0] - 2.0) < 1e-12, "symmetric_eigen eigenvalue");

  // A mixed pair: the library's own default on the left, TBB's on the right.
  lam::linalg::vector<double> plain({1.0, 1.0, 1.0});
  check(lam::linalg::dot(plain, v) == 7.0, "mixed dot with std::allocator");
  auto mixed = plain + v;
  static_assert(std::same_as<decltype(mixed), lam::linalg::vector<double>>);
  check(mixed[0] == 4.0, "mixed operator+ with a TBB allocator");

  std::println("{} container suite passed", name);
}

// cache_aligned_allocator's whole point: storage begins on a cache line.
void alignment_suite()
{
  using alloc = tbb::cache_aligned_allocator<double>;
  lam::linalg::vector<double, alloc> v(64, alloc{});

  auto addr = reinterpret_cast<std::uintptr_t>(v.as_span().data());
  // TBB guarantees cache-line alignment; 64 is the line size on every platform it supports here.
  check(addr % 64 == 0, "cache_aligned_allocator did not return cache-line-aligned storage");

  matrix<double, storage_layout::row_major, alloc> m(16, 16, alloc{});
  auto maddr = reinterpret_cast<std::uintptr_t>(&m[0, 0]);
  check(maddr % 64 == 0, "matrix storage not cache-line-aligned");

  std::println("alignment suite passed (vector at {:#x})", addr);
}

// TBB's real std::pmr resources need no adaptor; its feature macro wants libstdc++, so libc++ lacks them.
#if __TBB_CPP17_MEMORY_RESOURCE_PRESENT
void tbb_resource_suites()
{
  tbb::cache_aligned_resource aligned(std::pmr::new_delete_resource());

  failures += run_resource_suite("std::pmr::new_delete_resource (control)", std::pmr::new_delete_resource());
  failures += run_resource_suite("tbb::scalable_memory_resource", tbb::scalable_memory_resource());
  failures += run_resource_suite("tbb::cache_aligned_resource", &aligned);
}
#else
void tbb_resource_suites()
{ std::println("tbb pmr resources: unavailable in this build (TBB gates them on libstdc++)"); }
#endif

int main()
{
  std::println("allocator shapes:");
  allocator_shape<tbb::scalable_allocator<double>>("tbb::scalable_allocator");
  allocator_shape<tbb::cache_aligned_allocator<double>>("tbb::cache_aligned_allocator");

  container_suite<tbb::scalable_allocator<double>>("tbb::scalable_allocator");
  container_suite<tbb::cache_aligned_allocator<double>>("tbb::cache_aligned_allocator");
  alignment_suite();
  tbb_resource_suites();

  if (failures != 0)
  {
    std::println("{} failure(s)", failures);
    return 1;
  }
  std::println("\n✅ TBB allocator tests passed");
  return 0;
}
