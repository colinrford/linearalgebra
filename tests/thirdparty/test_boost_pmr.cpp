/*
 *  test_boost_pmr.cpp - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  test_boost_pmr is a c++ module
 *
 *  Correctness only: does lam.linalg behave when handed memory that came from
 *  outside the standard library? Boost has its own memory_resource hierarchy,
 *  unrelated by inheritance to std::pmr's, so reaching it at all takes an
 *  adaptor. Boost is optional -- this directory is skipped when it is absent.
 */

#include <boost/container/pmr/global_resource.hpp>
#include <boost/container/pmr/monotonic_buffer_resource.hpp>
#include <boost/container/pmr/unsynchronized_pool_resource.hpp>
#include <boost/pool/pool_alloc.hpp>

import std;
import lam.linearalgebra;

import linalg.test.pmr_suite;

using lam::linalg::matrix;
using lam::linalg::storage_layout;

int failures = 0;

// Prints only when nothing failed, so a summary cannot claim a pass over failed checks.
struct suite_guard
{
  int before = failures;
  std::string_view name;
  ~suite_guard()
  {
    if (failures == before)
      std::println("{} passed", name);
  }
};

void check(bool ok, std::string_view what)
{
  if (!ok)
  {
    std::println("FAILED: {}", what);
    ++failures;
  }
}

// The two memory_resource bases have the same shape, so the way across is to own one and be the other.
class boost_resource_adaptor final : public std::pmr::memory_resource
{
public:
  explicit boost_resource_adaptor(boost::container::pmr::memory_resource& upstream) : m_upstream{&upstream} {}

  std::size_t allocations() const { return m_allocations; }

private:
  void* do_allocate(std::size_t bytes, std::size_t alignment) override
  {
    ++m_allocations;
    return m_upstream->allocate(bytes, alignment);
  }

  void do_deallocate(void* p, std::size_t bytes, std::size_t alignment) override
  { m_upstream->deallocate(p, bytes, alignment); }

  bool do_is_equal(const std::pmr::memory_resource& other) const noexcept override
  {
    const auto* rhs = dynamic_cast<const boost_resource_adaptor*>(&other);
    return rhs != nullptr && rhs->m_upstream == m_upstream;
  }

  boost::container::pmr::memory_resource* m_upstream;
  std::size_t m_allocations = 0;
};

// Buffer-backed only: our storage must land in the caller's buffer, and upstream is the null resource.
template<typename Resource, typename... Upstream>
void buffer_residency_suite(std::string_view name, Upstream... upstream)
{
  suite_guard guard{.name = name};

  alignas(std::max_align_t) std::array<std::byte, 8192> buffer{};
  Resource arena(buffer.data(), buffer.size(), upstream...);
  boost_resource_adaptor adaptor(arena);

  const auto* lo = buffer.data();
  const auto* hi = buffer.data() + buffer.size();
  auto inside = [&](const void* p) {
    const auto* c = static_cast<const std::byte*>(p);
    return c >= lo && c < hi;
  };

  lam::linalg::pmr::vector<double> v({3.0, 4.0, 0.0}, &adaptor);
  check(inside(v.as_span().data()), "vector storage did not come from the buffer");

  lam::linalg::pmr::matrix<double> m(8, 8, &adaptor);
  for (std::size_t i = 0; i < 8; ++i)
    m[i, i] = 2.0;
  check(inside(&m[0, 0]), "matrix storage did not come from the buffer");

  // If result propagation is broken this address lands on the heap instead.
  auto product = m * m;
  check(inside(&product[0, 0]), "operator* result escaped the buffer");
  check(adaptor.allocations() >= 3, "adaptor saw fewer allocations than containers created");
}

// Boost.Pool's plain Allocator needs no adaptor, and is a singleton pool, so instances compare equal.
void boost_pool_allocator_suite()
{
  suite_guard guard{.name = "boost::pool_allocator"};

  using alloc = boost::pool_allocator<double>;
  using traits = std::allocator_traits<alloc>;
  static_assert(std::same_as<typename traits::value_type, double>);
  static_assert(std::same_as<typename traits::pointer, double*>, "lam.linalg stores raw T*, not a fancy pointer");
  static_assert(std::same_as<typename traits::template rebind_alloc<int>::value_type, int>);
  static_assert(traits::is_always_equal::value, "boost::pool_allocator is a singleton pool");

  lam::linalg::vector<double, alloc> v({3.0, 4.0, 0.0}, alloc{});
  check(v.norm() == 5.0, "pool_allocator vector norm");

  lam::linalg::vector<double, alloc> w = v + v;
  check(w[1] == 8.0, "pool_allocator operator+");
  check(w.get_allocator() == v.get_allocator(), "result kept the operand's allocator");

  matrix<double, storage_layout::row_major, alloc> m(4, 4, alloc{});
  for (std::size_t i = 0; i < 4; ++i)
    m[i, i] = 3.0;
  auto squared = m * m;
  check(squared[2, 2] == 9.0, "pool_allocator operator*");
}

int main()
{
  // The same suite against std and Boost resources; the std rows are the control.
  std::pmr::monotonic_buffer_resource std_arena;
  std::pmr::unsynchronized_pool_resource std_pool;
  boost::container::pmr::monotonic_buffer_resource boost_arena;
  boost::container::pmr::unsynchronized_pool_resource boost_pool;
  boost_resource_adaptor boost_arena_mr(boost_arena);
  boost_resource_adaptor boost_pool_mr(boost_pool);

  failures += run_resource_suite("std::pmr::new_delete_resource", std::pmr::new_delete_resource());
  failures += run_resource_suite("std::pmr::monotonic_buffer_resource", &std_arena);
  failures += run_resource_suite("std::pmr::unsynchronized_pool_resource", &std_pool);
  failures += run_resource_suite("boost monotonic_buffer_resource", &boost_arena_mr);
  failures += run_resource_suite("boost unsynchronized_pool_resource", &boost_pool_mr);

  buffer_residency_suite<boost::container::pmr::monotonic_buffer_resource>(
    "boost buffer residency", boost::container::pmr::null_memory_resource());

  boost_pool_allocator_suite();

  if (failures != 0)
  {
    std::println("{} failure(s)", failures);
    return 1;
  }
  std::println("\n✅ Boost pmr/allocator tests passed");
  return 0;
}
