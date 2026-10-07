/*
 *  test_allocator_traits.cpp - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  test_allocator_traits is a c++ module
 */

import std;
import lam.linearalgebra;

using lam::linalg::vector;

int failures = 0;

void check(bool ok, std::string_view what)
{
  if (!ok)
  {
    std::println("FAILED: {}", what);
    ++failures;
  }
}

template<typename T, bool Pocca, bool Pocs, bool Soccc>
struct traits_alloc
{
  using value_type = T;
  using propagate_on_container_copy_assignment = std::bool_constant<Pocca>;
  using propagate_on_container_swap = std::bool_constant<Pocs>;
  using is_always_equal = std::false_type;
  int tag = 0;

  constexpr traits_alloc() = default;
  constexpr explicit traits_alloc(int t) : tag{t} {}
  template<typename U>
  constexpr traits_alloc(const traits_alloc<U, Pocca, Pocs, Soccc>& o) : tag{o.tag}
  {}

  constexpr T* allocate(std::size_t n) { return std::allocator<T>{}.allocate(n); }
  constexpr void deallocate(T* p, std::size_t n) { std::allocator<T>{}.deallocate(p, n); }

  constexpr traits_alloc select_on_container_copy_construction() const
    requires Soccc
  { return traits_alloc{0}; }

  constexpr bool operator==(const traits_alloc& o) const { return tag == o.tag; }
};

template<typename T>
using pocca_alloc = traits_alloc<T, true, false, false>;
template<typename T>
using soccc_alloc = traits_alloc<T, false, false, true>;

// std::allocator propagates on move only; pmr propagates on nothing and customizes SOCCC.
void trait_shapes()
{
  using std_traits = std::allocator_traits<std::allocator<double>>;
  static_assert(std_traits::propagate_on_container_move_assignment::value);
  static_assert(!std_traits::propagate_on_container_copy_assignment::value);
  static_assert(!std_traits::propagate_on_container_swap::value);

  using pmr_traits = std::allocator_traits<std::pmr::polymorphic_allocator<double>>;
  static_assert(!pmr_traits::propagate_on_container_move_assignment::value);
  static_assert(!pmr_traits::propagate_on_container_copy_assignment::value);
  static_assert(!pmr_traits::propagate_on_container_swap::value);
}

// lam keeps the target's allocator on copy-assignment even when POCCA asks for the source's.
void copy_assignment_ignores_pocca()
{
  vector<double, pocca_alloc<double>> target({1.0, 2.0}, pocca_alloc<double>{1});
  const vector<double, pocca_alloc<double>> source({3.0, 4.0}, pocca_alloc<double>{2});
  target = source;
  check(target[0] == 3.0 && target[1] == 4.0, "copy-assign values");
  check(target.get_allocator().tag == 1, "copy-assign keeps the target's allocator under POCCA");
}

// lam copies the source's allocator; std would call select_on_container_copy_construction.
void copy_construction_ignores_soccc()
{
  const vector<double, soccc_alloc<double>> source({1.0, 2.0}, soccc_alloc<double>{7});
  const auto copy = source;
  check(copy.get_allocator().tag == 7, "copy-construct keeps the source's allocator, not SOCCC's");
}

// Records which allocator made each block, so a swap that frees through the wrong one is caught.
std::map<const void*, int> block_owner;
int allocations = 0;
int wrong_frees = 0;

template<typename T, bool Pocs>
struct owned_alloc
{
  using value_type = T;
  using propagate_on_container_swap = std::bool_constant<Pocs>;
  using is_always_equal = std::false_type;
  template<typename U>
  struct rebind
  {
    using other = owned_alloc<U, Pocs>;
  };
  int tag = 0;

  owned_alloc() = default;
  explicit owned_alloc(int t) : tag{t} {}
  template<typename U>
  owned_alloc(const owned_alloc<U, Pocs>& o) : tag{o.tag}
  {}

  T* allocate(std::size_t n)
  {
    ++allocations;
    T* p = std::allocator<T>{}.allocate(n);
    block_owner[p] = tag;
    return p;
  }
  void deallocate(T* p, std::size_t n)
  {
    if (block_owner[p] != tag)
      ++wrong_frees;
    block_owner.erase(p);
    std::allocator<T>{}.deallocate(p, n);
  }

  bool operator==(const owned_alloc& o) const { return tag == o.tag; }
};

// POCS: allocators travel with their storage. Otherwise swap refuses and leaves both untouched.
template<template<typename, typename> class V, bool Pocs>
void swap_unequal_allocators(std::string_view name)
{
  using A = owned_alloc<double, Pocs>;
  wrong_frees = 0;
  {
    V<double, A> a({1.0, 2.0}, A{1});
    V<double, A> b({3.0, 4.0, 5.0}, A{2});
    allocations = 0;
    using std::swap;
    if constexpr (Pocs)
    {
      swap(a, b);
      check(a.size() == 3 && a[2] == 5.0 && b.size() == 2 && b[1] == 2.0, std::format("{}: swap values", name));
      check(a.get_allocator().tag == 2 && b.get_allocator().tag == 1, std::format("{}: POCS swaps allocators", name));
    }
    else
    {
      bool threw = false;
      try
      {
        swap(a, b);
      }
      catch (const lam::linalg::vector_exception::unequal_allocators&)
      {
        threw = true;
      }
      check(threw, std::format("{}: non-POCS swap of unequal allocators throws", name));
      check(a.size() == 2 && a[1] == 2.0 && b.size() == 3 && b[2] == 5.0, std::format("{}: refused swap", name));
      check(a.get_allocator().tag == 1 && b.get_allocator().tag == 2, std::format("{}: refused swap", name));
    }
    check(allocations == 0, std::format("{}: swap does not allocate", name));
  }
  check(wrong_frees == 0, std::format("{}: every block freed by the allocator that made it", name));
}

template<template<typename, typename> class V>
void swap_equal_allocators(std::string_view name)
{
  using A = owned_alloc<double, false>;
  V<double, A> a({1.0, 2.0}, A{1});
  V<double, A> b({3.0, 4.0, 5.0}, A{1});
  allocations = 0;
  std::ranges::swap(a, b);
  check(a.size() == 3 && b.size() == 2 && a[0] == 3.0 && b[0] == 1.0, std::format("{}: equal swap values", name));
  check(allocations == 0, std::format("{}: equal-allocator swap does not allocate", name));
  swap(a, a);
  check(a.size() == 3 && a[1] == 4.0, std::format("{}: self-swap is a no-op", name));
}

// noexcept exactly when swap can always exchange pointers.
template<template<typename, typename> class V>
consteval bool swap_noexcept_shapes()
{
  static_assert(std::is_nothrow_swappable_v<V<double, std::allocator<double>>>);
  static_assert(std::is_nothrow_swappable_v<V<double, owned_alloc<double, true>>>);
  static_assert(!std::is_nothrow_swappable_v<V<double, owned_alloc<double, false>>>);
  static_assert(!std::is_nothrow_swappable_v<V<double, std::pmr::polymorphic_allocator<double>>>);
  return true;
}

// POCS with unequal allocators, or non-POCS with equal ones; the refused case would not be a constant expression.
template<template<typename, typename> class V, bool Pocs>
consteval bool constexpr_swap()
{
  using A = traits_alloc<double, false, Pocs, false>;
  const int b_tag = Pocs ? 2 : 1;
  V<double, A> a({1.0, 2.0}, A{1});
  V<double, A> b({3.0, 4.0, 5.0}, A{b_tag});
  using std::swap;
  swap(a, b);
  return a.size() == 3 && a[2] == 5.0 && b[1] == 2.0 && a.get_allocator().tag == b_tag && b.get_allocator().tag == 1;
}

template<typename A>
using row_mat = lam::linalg::matrix<double, lam::linalg::storage_layout::row_major, A>;

consteval bool matrix_swap_noexcept_shapes()
{
  static_assert(std::is_nothrow_swappable_v<row_mat<std::allocator<double>>>);
  static_assert(std::is_nothrow_swappable_v<row_mat<owned_alloc<double, true>>>);
  static_assert(!std::is_nothrow_swappable_v<row_mat<owned_alloc<double, false>>>);
  static_assert(!std::is_nothrow_swappable_v<row_mat<std::pmr::polymorphic_allocator<double>>>);
  return true;
}

template<bool Pocs>
consteval bool constexpr_matrix_swap()
{
  using A = traits_alloc<double, false, Pocs, false>;
  const int b_tag = Pocs ? 2 : 1;
  row_mat<A> a({{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}}, A{1});
  row_mat<A> b({{7.0, 8.0}, {9.0, 10.0}, {11.0, 12.0}}, A{b_tag});
  using std::swap;
  swap(a, b);
  return a.rows() == 3 && a.cols() == 2 && a[2, 1] == 12.0 && b.rows() == 2 && b.cols() == 3 && b[1, 2] == 6.0 &&
         a.get_allocator().tag == b_tag && b.get_allocator().tag == 1;
}

static_assert(swap_noexcept_shapes<lam::linalg::vector>());
static_assert(swap_noexcept_shapes<lam::linalg::r_vector>());
static_assert(swap_noexcept_shapes<lam::linalg::raw_vector>());
static_assert(matrix_swap_noexcept_shapes());
static_assert(constexpr_swap<lam::linalg::vector, true>());
static_assert(constexpr_swap<lam::linalg::vector, false>());
static_assert(constexpr_swap<lam::linalg::r_vector, true>());
static_assert(constexpr_swap<lam::linalg::r_vector, false>());
static_assert(constexpr_swap<lam::linalg::raw_vector, true>());
static_assert(constexpr_swap<lam::linalg::raw_vector, false>());
static_assert(constexpr_matrix_swap<true>());
static_assert(constexpr_matrix_swap<false>());

template<bool Pocs>
void matrix_swap(std::string_view name)
{
  using A = owned_alloc<double, Pocs>;
  using mat = lam::linalg::matrix<double, lam::linalg::storage_layout::row_major, A>;
  wrong_frees = 0;
  {
    mat a({{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}}, A{1});
    mat b({{7.0, 8.0}, {9.0, 10.0}, {11.0, 12.0}}, A{2});
    allocations = 0;
    using std::swap;
    if constexpr (Pocs)
    {
      swap(a, b);
      check(a.rows() == 3 && a.cols() == 2 && a[2, 1] == 12.0, std::format("{}: swapped shape and values (a)", name));
      check(b.rows() == 2 && b.cols() == 3 && b[1, 2] == 6.0, std::format("{}: swapped shape and values (b)", name));
      check(a.get_allocator().tag == 2 && b.get_allocator().tag == 1, std::format("{}: POCS swaps allocators", name));
    }
    else
    {
      bool threw = false;
      try
      {
        swap(a, b);
      }
      catch (const lam::linalg::matrix_exception::unequal_allocators&)
      {
        threw = true;
      }
      check(threw, std::format("{}: non-POCS swap of unequal allocators throws", name));
      check(a.rows() == 2 && a.cols() == 3 && a[1, 2] == 6.0, std::format("{}: refused swap (a)", name));
      check(b.rows() == 3 && b.cols() == 2 && b[2, 1] == 12.0, std::format("{}: refused swap (b)", name));
      check(a.get_allocator().tag == 1 && b.get_allocator().tag == 2, std::format("{}: refused swap", name));
    }
    check(allocations == 0, std::format("{}: swap does not allocate", name));
  }
  check(wrong_frees == 0, std::format("{}: every block freed by the allocator that made it", name));
}

int main()
{
  trait_shapes();
  copy_assignment_ignores_pocca();
  copy_construction_ignores_soccc();

  swap_unequal_allocators<lam::linalg::vector, true>("vector POCS");
  swap_unequal_allocators<lam::linalg::vector, false>("vector non-POCS");
  swap_unequal_allocators<lam::linalg::r_vector, true>("r_vector POCS");
  swap_unequal_allocators<lam::linalg::r_vector, false>("r_vector non-POCS");
  swap_unequal_allocators<lam::linalg::raw_vector, true>("raw_vector POCS");
  swap_unequal_allocators<lam::linalg::raw_vector, false>("raw_vector non-POCS");
  swap_equal_allocators<lam::linalg::vector>("vector");
  swap_equal_allocators<lam::linalg::r_vector>("r_vector");
  swap_equal_allocators<lam::linalg::raw_vector>("raw_vector");
  matrix_swap<true>("matrix POCS");
  matrix_swap<false>("matrix non-POCS");

  if (failures != 0)
  {
    std::println("{} failure(s)", failures);
    return 1;
  }
  std::println("\n✅ allocator trait tests passed");
  return 0;
}
