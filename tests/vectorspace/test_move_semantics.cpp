/*
 *  test_move_semantics.cpp - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  test_move_semantics is a c++ module
 */

import std;
import lam.linearalgebra;

using lam::linalg::storage_layout;

// Stateful allocator usable in constant evaluation; each instance counts its own live blocks.
template<typename T, bool Pocma>
struct tagged_alloc
{
  using value_type = T;
  using propagate_on_container_move_assignment = std::bool_constant<Pocma>;
  template<typename U>
  struct rebind
  {
    using other = tagged_alloc<U, Pocma>;
  };

  int id;
  int* live;

  constexpr tagged_alloc(int i, int* l) noexcept : id{i}, live{l} {}
  template<typename U>
  constexpr tagged_alloc(const tagged_alloc<U, Pocma>& o) noexcept : id{o.id}, live{o.live}
  {}

  constexpr T* allocate(std::size_t n)
  {
    ++*live;
    return std::allocator<T>{}.allocate(n);
  }
  constexpr void deallocate(T* p, std::size_t n)
  {
    --*live;
    std::allocator<T>{}.deallocate(p, n);
  }

  template<typename U>
  friend constexpr bool operator==(const tagged_alloc& a, const tagged_alloc<U, Pocma>& b) noexcept
  { return a.id == b.id; }
};

// Stateless C++11-style allocator: always equal, does not propagate. Counts live blocks at runtime.
template<typename T>
struct stateless_alloc
{
  using value_type = T;
  using propagate_on_container_move_assignment = std::false_type;
  using is_always_equal = std::true_type;

  static inline int live = 0;

  constexpr stateless_alloc() noexcept = default;
  template<typename U>
  constexpr stateless_alloc(const stateless_alloc<U>&) noexcept
  {}

  constexpr T* allocate(std::size_t n)
  {
    if !consteval
    {
      ++live;
    }
    return std::allocator<T>{}.allocate(n);
  }
  constexpr void deallocate(T* p, std::size_t n)
  {
    if !consteval
    {
      --live;
    }
    std::allocator<T>{}.deallocate(p, n);
  }

  template<typename U>
  friend constexpr bool operator==(const stateless_alloc&, const stateless_alloc<U>&) noexcept
  { return true; }
};

// Counts live allocations; throws bad_alloc once `allow` runs out.
struct counting_resource : std::pmr::memory_resource
{
  int live = 0;
  std::size_t allow = std::numeric_limits<std::size_t>::max();
  void* do_allocate(std::size_t n, std::size_t al) override
  {
    if (allow == 0)
      throw std::bad_alloc();
    --allow;
    ++live;
    return std::pmr::new_delete_resource()->allocate(n, al);
  }
  void do_deallocate(void* p, std::size_t n, std::size_t al) override
  {
    --live;
    std::pmr::new_delete_resource()->deallocate(p, n, al);
  }
  bool do_is_equal(const std::pmr::memory_resource& o) const noexcept override { return this == &o; }
};

using pmr_alloc = std::pmr::polymorphic_allocator<double>;

// Returns 0 on success, else the number of the first failed check. a1 and a2 may or may not compare equal.
template<template<typename, typename> class V, typename A>
constexpr int vector_move_suite(const A& a1, const A& a2)
{
  using vec = V<double, A>;
  const bool pocma = std::allocator_traits<A>::propagate_on_container_move_assignment::value;
  const bool steals = pocma || a1 == a2;

  vec src({1.0, 2.0, 3.0}, a1);
  vec moved(std::move(src));
  if (src.size() != 0)
    return 1;
  if (moved.size() != 3 || moved[2] != 3.0 || !(moved.get_allocator() == a1))
    return 2;

  vec from_empty({4.0, 5.0}, a2);
  from_empty = std::move(src);
  if (from_empty.size() != 0 || !(from_empty.get_allocator() == (pocma ? a1 : a2)))
    return 3;

  vec equal_target({1.0, 2.0}, a1);
  vec equal_source({7.0, 8.0, 9.0}, a1);
  equal_target = std::move(equal_source);
  if (equal_target.size() != 3 || equal_target[0] != 7.0 || equal_source.size() != 0)
    return 4;

  vec same_size({1.0, 2.0}, a2);
  vec same_size_source({5.0, 6.0}, a1);
  same_size = std::move(same_size_source);
  if (same_size.size() != 2 || same_size[1] != 6.0 || !(same_size.get_allocator() == (pocma ? a1 : a2)))
    return 5;
  if (same_size_source.size() != (steals ? 0 : 2))
    return 6;

  vec grow({1.0, 2.0}, a2);
  vec grow_source({5.0, 6.0, 7.0}, a1);
  grow = std::move(grow_source);
  if (grow.size() != 3 || grow[2] != 7.0 || !(grow.get_allocator() == (pocma ? a1 : a2)))
    return 7;
  if (grow_source.size() != (steals ? 0 : 3))
    return 8;

  vec copy_target({1.0, 2.0}, a2);
  copy_target = moved;
  if (copy_target.size() != 3 || copy_target[1] != 2.0 || !(copy_target.get_allocator() == a2) || moved.size() != 3)
    return 9;

  vec copied(moved);
  if (copied.size() != 3 || copied[0] != 1.0 || !(copied.get_allocator() == a1))
    return 10;

  vec& alias = moved;
  moved = std::move(alias);
  if (moved.size() != 3 || moved[0] != 1.0)
    return 11;

  src = vec({8.0, 9.0}, a1);
  if (src.size() != 2 || src[1] != 9.0)
    return 12;

  vec z = zero_like(src);
  if (z.size() != src.size() || !(z.get_allocator() == src.get_allocator()) || z[0] != 0.0)
    return 14;

  vec twice(std::move(moved));
  vec after_twice({1.0, 2.0}, a1);
  after_twice = std::move(moved);
  if (after_twice.size() != 0 || twice.size() != 3)
    return 13;

  return 0;
}

template<typename A>
constexpr int matrix_move_suite(const A& a1, const A& a2)
{
  using mat = lam::linalg::matrix<double, storage_layout::row_major, A>;
  const bool pocma = std::allocator_traits<A>::propagate_on_container_move_assignment::value;
  const bool steals = pocma || a1 == a2;

  mat src(2, 3, a1);
  src[1, 2] = 4.0;
  mat moved(std::move(src));
  if (src.rows() != 0 || src.cols() != 0)
    return 1;
  if (moved.rows() != 2 || moved.cols() != 3 || moved[1, 2] != 4.0)
    return 2;

  mat from_empty(2, 2, a2);
  from_empty = std::move(src);
  if (from_empty.rows() != 0 || from_empty.cols() != 0)
    return 3;

  mat same_size(2, 3, a2);
  same_size = std::move(moved);
  if (same_size[1, 2] != 4.0 || !(same_size.get_allocator() == (pocma ? a1 : a2)))
    return 4;
  const std::size_t expected_rows = steals ? 0 : 2;
  const std::size_t expected_cols = steals ? 0 : 3;
  if (moved.rows() != expected_rows || moved.cols() != expected_cols)
    return 5;

  mat reshape(3, 1, a2);
  reshape = std::move(same_size);
  if (reshape.rows() != 2 || reshape.cols() != 3 || reshape[1, 2] != 4.0)
    return 6;

  mat copy_target(1, 1, a2);
  copy_target = reshape;
  if (copy_target.rows() != 2 || copy_target[1, 2] != 4.0 || !(copy_target.get_allocator() == a2))
    return 7;

  mat& alias = reshape;
  reshape = std::move(alias);
  if (reshape.rows() != 2 || reshape[1, 2] != 4.0)
    return 8;

  auto transposed = reshape.transpose();
  if (transposed.rows() != 3 || transposed[2, 1] != 4.0 || !(transposed.get_allocator() == reshape.get_allocator()))
    return 9;

  auto negated = -reshape;
  if (negated[1, 2] != -4.0 || !(negated.get_allocator() == reshape.get_allocator()))
    return 10;

  auto mz = zero_like(reshape);
  if (mz.rows() != reshape.rows() || mz.cols() != reshape.cols() || mz[1, 2] != 0.0 ||
      !(mz.get_allocator() == reshape.get_allocator()))
    return 12;

  auto mi = identity_like(reshape.transpose() * reshape);
  if (mi.rows() != 3 || mi[0, 0] != 1.0 || mi[0, 1] != 0.0 || !(mi.get_allocator() == reshape.get_allocator()))
    return 13;

  using float_alloc = typename std::allocator_traits<A>::template rebind_alloc<float>;
  if (!(A(float_alloc(a1)) == a1))
    return 14;

  return 0;
}

// On allocation failure, assignment leaves both sides untouched.
template<template<typename, typename> class V>
int vector_strong_guarantee()
{
  counting_resource r1;
  counting_resource r2;
  {
    using vec = V<double, pmr_alloc>;
    vec target({1.0, 2.0}, pmr_alloc{&r2});
    vec source({5.0, 6.0, 7.0}, pmr_alloc{&r1});
    r2.allow = 0;
    try
    {
      target = std::move(source);
      return 1;
    }
    catch (const std::bad_alloc&)
    {}
    if (target.size() != 2 || target[1] != 2.0 || source.size() != 3 || source[2] != 7.0)
      return 2;
    try
    {
      target = source;
      return 3;
    }
    catch (const std::bad_alloc&)
    {}
    if (target.size() != 2 || target[0] != 1.0 || source.size() != 3)
      return 4;
    r2.allow = std::numeric_limits<std::size_t>::max();
    target = std::move(source);
    if (target.size() != 3 || target[2] != 7.0)
      return 5;
  }
  return (r1.live || r2.live) ? 100 : 0;
}

int matrix_strong_guarantee()
{
  counting_resource r1;
  counting_resource r2;
  {
    using mat = lam::linalg::matrix<double, storage_layout::row_major, pmr_alloc>;
    mat target(2, 2, pmr_alloc{&r2});
    target[1, 1] = 3.0;
    mat source(2, 3, pmr_alloc{&r1});
    r2.allow = 0;
    try
    {
      target = std::move(source);
      return 1;
    }
    catch (const std::bad_alloc&)
    {}
    if (target.rows() != 2 || target.cols() != 2 || target[1, 1] != 3.0 || source.rows() != 2 || source.cols() != 3)
      return 2;
  }
  return (r1.live || r2.live) ? 100 : 0;
}

template<template<typename, typename> class V, bool Pocma>
constexpr int tagged_vector_run()
{
  int live1 = 0;
  int live2 = 0;
  int result = vector_move_suite<V>(tagged_alloc<double, Pocma>{1, &live1}, tagged_alloc<double, Pocma>{2, &live2});
  if (result != 0)
    return result;
  if (live1 != 0 || live2 != 0)
    return 100;
  return 0;
}

template<bool Pocma>
constexpr int tagged_matrix_run()
{
  int live1 = 0;
  int live2 = 0;
  int result = matrix_move_suite(tagged_alloc<double, Pocma>{1, &live1}, tagged_alloc<double, Pocma>{2, &live2});
  if (result != 0)
    return result;
  if (live1 != 0 || live2 != 0)
    return 100;
  return 0;
}

template<template<typename, typename> class V>
constexpr int std_vector_run()
{ return vector_move_suite<V>(std::allocator<double>{}, std::allocator<double>{}); }

constexpr int std_matrix_run() { return matrix_move_suite(std::allocator<double>{}, std::allocator<double>{}); }

template<template<typename, typename> class V>
constexpr int stateless_vector_run()
{
  if !consteval
  {
    stateless_alloc<double>::live = 0;
  }
  int result = vector_move_suite<V>(stateless_alloc<double>{}, stateless_alloc<double>{});
  if !consteval
  {
    if (result == 0 && stateless_alloc<double>::live != 0)
      return 100;
  }
  return result;
}

constexpr int stateless_matrix_run()
{
  if !consteval
  {
    stateless_alloc<double>::live = 0;
  }
  int result = matrix_move_suite(stateless_alloc<double>{}, stateless_alloc<double>{});
  if !consteval
  {
    if (result == 0 && stateless_alloc<double>::live != 0)
      return 100;
  }
  return result;
}

static_assert(tagged_vector_run<lam::linalg::vector, false>() == 0);
static_assert(tagged_vector_run<lam::linalg::vector, true>() == 0);
static_assert(std_vector_run<lam::linalg::vector>() == 0);
static_assert(stateless_vector_run<lam::linalg::vector>() == 0);
static_assert(tagged_vector_run<lam::linalg::r_vector, false>() == 0);
static_assert(tagged_vector_run<lam::linalg::r_vector, true>() == 0);
static_assert(std_vector_run<lam::linalg::r_vector>() == 0);
static_assert(stateless_vector_run<lam::linalg::r_vector>() == 0);
static_assert(tagged_vector_run<lam::linalg::raw_vector, false>() == 0);
static_assert(tagged_vector_run<lam::linalg::raw_vector, true>() == 0);
static_assert(std_vector_run<lam::linalg::raw_vector>() == 0);
static_assert(stateless_vector_run<lam::linalg::raw_vector>() == 0);
static_assert(tagged_matrix_run<false>() == 0);
static_assert(tagged_matrix_run<true>() == 0);
static_assert(std_matrix_run() == 0);
static_assert(stateless_matrix_run() == 0);

// Runs the suites over two allocators built from `make(resource)`; checks every resource's upstream is released.
template<typename Resource, typename Make>
int pmr_resource_run(std::string_view what, Make make)
{
  counting_resource up1;
  counting_resource up2;
  int result = 0;
  {
    Resource res1{&up1};
    Resource res2{&up2};
    auto a1 = make(res1);
    auto a2 = make(res2);

    // Each suite runs twice: with two different resources, then with one shared resource.
    if ((result = vector_move_suite<lam::linalg::vector>(a1, a2)) == 0)
      if ((result = vector_move_suite<lam::linalg::vector>(a1, a1)) == 0)
        if ((result = vector_move_suite<lam::linalg::r_vector>(a1, a2)) == 0)
          if ((result = vector_move_suite<lam::linalg::r_vector>(a1, a1)) == 0)
            if ((result = vector_move_suite<lam::linalg::raw_vector>(a1, a2)) == 0)
              if ((result = vector_move_suite<lam::linalg::raw_vector>(a1, a1)) == 0)
                if ((result = matrix_move_suite(a1, a2)) == 0)
                  result = matrix_move_suite(a1, a1);
  }
  if (result == 0 && (up1.live || up2.live))
    result = 100;
  if (result != 0)
    std::println("{}: check {} failed", what, result);
  return result;
}

struct counting_as_upstream : counting_resource
{
  explicit counting_as_upstream(counting_resource* up) : upstream{up} {}
  counting_resource* upstream;
  void* do_allocate(std::size_t n, std::size_t al) override { return upstream->allocate(n, al); }
  void do_deallocate(void* p, std::size_t n, std::size_t al) override { upstream->deallocate(p, n, al); }
};

int failures = 0;

void check(int result, std::string_view what)
{
  if (result != 0)
  {
    std::println("FAILED {}: check {}", what, result);
    ++failures;
  }
  else
    std::println("{} passed", what);
}

int main()
{
  using lam::linalg::r_vector;
  using lam::linalg::raw_vector;
  using lam::linalg::vector;

  check(tagged_vector_run<vector, false>(), "vector, stateful non-propagating");
  check(tagged_vector_run<vector, true>(), "vector, stateful propagating");
  check(std_vector_run<vector>(), "vector, std::allocator");
  check(stateless_vector_run<vector>(), "vector, stateless always-equal non-propagating");
  check(tagged_vector_run<r_vector, false>(), "r_vector, stateful non-propagating");
  check(tagged_vector_run<r_vector, true>(), "r_vector, stateful propagating");
  check(std_vector_run<r_vector>(), "r_vector, std::allocator");
  check(stateless_vector_run<r_vector>(), "r_vector, stateless always-equal non-propagating");
  check(tagged_vector_run<raw_vector, false>(), "raw_vector, stateful non-propagating");
  check(tagged_vector_run<raw_vector, true>(), "raw_vector, stateful propagating");
  check(std_vector_run<raw_vector>(), "raw_vector, std::allocator");
  check(stateless_vector_run<raw_vector>(), "raw_vector, stateless always-equal non-propagating");
  check(tagged_matrix_run<false>(), "matrix, stateful non-propagating");
  check(tagged_matrix_run<true>(), "matrix, stateful propagating");
  check(std_matrix_run(), "matrix, std::allocator");
  check(stateless_matrix_run(), "matrix, stateless always-equal non-propagating");

  check(pmr_resource_run<counting_as_upstream>("pmr, plain resource", [](auto& r) { return pmr_alloc{&r}; }),
        "all types, pmr plain resource (equal and unequal)");
  check(pmr_resource_run<std::pmr::monotonic_buffer_resource>("pmr, monotonic", [](auto& r) { return pmr_alloc{&r}; }),
        "all types, pmr monotonic_buffer_resource");
  check(pmr_resource_run<std::pmr::unsynchronized_pool_resource>("pmr, unsynchronized pool",
                                                                 [](auto& r) { return pmr_alloc{&r}; }),
        "all types, pmr unsynchronized_pool_resource");
  check(pmr_resource_run<std::pmr::synchronized_pool_resource>("pmr, synchronized pool",
                                                               [](auto& r) { return pmr_alloc{&r}; }),
        "all types, pmr synchronized_pool_resource");
  check(pmr_resource_run<counting_as_upstream>(
          "scoped_allocator_adaptor", [](auto& r) { return std::scoped_allocator_adaptor<pmr_alloc>{pmr_alloc{&r}}; }),
        "all types, scoped_allocator_adaptor<polymorphic_allocator>");

  check(vector_strong_guarantee<vector>(), "vector, strong guarantee on allocation failure");
  check(vector_strong_guarantee<r_vector>(), "r_vector, strong guarantee on allocation failure");
  check(vector_strong_guarantee<raw_vector>(), "raw_vector, strong guarantee on allocation failure");
  check(matrix_strong_guarantee(), "matrix, strong guarantee on allocation failure");

  if (failures != 0)
  {
    std::println("{} failure(s)", failures);
    return 1;
  }
  return 0;
}
