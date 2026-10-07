/*
 *  test_mixed_allocators.cpp - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  test_mixed_allocators is a c++ module
 */

import std;
import lam.linearalgebra;

using lam::linalg::storage_layout;

// Stateful and distinct from std::allocator and polymorphic_allocator, so the rule is observable.
template<typename T>
struct tagged_alloc
{
  using value_type = T;
  int tag = 0;

  constexpr tagged_alloc() = default;
  constexpr explicit tagged_alloc(int t) : tag{t} {}
  template<typename U>
  constexpr tagged_alloc(const tagged_alloc<U>& o) : tag{o.tag}
  {}

  constexpr T* allocate(std::size_t n) { return std::allocator<T>{}.allocate(n); }
  constexpr void deallocate(T* p, std::size_t n) { std::allocator<T>{}.deallocate(p, n); }

  constexpr bool operator==(const tagged_alloc& o) const { return tag == o.tag; }
};

int failures = 0;

void check(bool ok, std::string_view what)
{
  if (!ok)
  {
    std::println("FAILED: {}", what);
    ++failures;
  }
}

template<template<typename, typename> class V>
void mixed_suite(std::string_view name)
{
  using std_v = V<double, std::allocator<double>>;
  using tag_v = V<double, tagged_alloc<double>>;
  using pmr_v = V<double, std::pmr::polymorphic_allocator<double>>;

  std::array<std::byte, 8192> buffer;
  std::pmr::monotonic_buffer_resource pool{buffer.data(), buffer.size()};

  const std_v a({1.0, 2.0, 3.0});
  const tag_v b({10.0, 20.0, 30.0}, tagged_alloc<double>{7});
  const pmr_v c({100.0, 200.0, 300.0}, &pool);

  // The result type is the LEFT operand's, allocator included.
  static_assert(std::same_as<decltype(a + b), std_v>);
  static_assert(std::same_as<decltype(b + a), tag_v>);
  static_assert(std::same_as<decltype(c + a), pmr_v>);

  {
    auto s = a + b;
    check(s[0] == 11.0 && s[2] == 33.0, "mixed operator+ values");
  }
  {
    auto s = b + a;
    check(s[0] == 11.0 && s[2] == 33.0, "mixed operator+ values (reversed)");
    check(s.get_allocator().tag == 7, "mixed operator+ must keep the left operand's allocator");
  }
  {
    auto d = c - a;
    check(d[0] == 99.0 && d[2] == 297.0, "mixed operator- values");
    check(d.get_allocator().resource() == &pool, "mixed operator- must keep the left operand's resource");
  }

  // Compound assignment leaves the target's allocator alone, by construction.
  {
    std_v t({1.0, 1.0, 1.0});
    t += b;
    check(t[1] == 21.0, "mixed operator+= values");
    t -= c;
    check(t[0] == 11.0 - 100.0, "mixed operator-= values");
  }

  // Scalar-valued algorithms do not care about allocators at all.
  check(lam::linalg::dot(a, b) == 1.0 * 10.0 + 2.0 * 20.0 + 3.0 * 30.0, "mixed dot");
  check(a.dot(c) == 100.0 + 400.0 + 900.0, "mixed member dot");
  check(lam::linalg::distance(a, a) == 0.0, "mixed distance with itself");
  check(a == std_v({1.0, 2.0, 3.0}), "same-type equality still works");
  check(!(a == b), "mixed equality compares values");
  {
    const tag_v same({1.0, 2.0, 3.0}, tagged_alloc<double>{9});
    check(a == same, "mixed equality is true for equal values on different allocators");

    // Different lengths are unequal, and comparing them must not read past the shorter one.
    const tag_v shorter({1.0, 2.0}, tagged_alloc<double>{9});
    const tag_v longer({1.0, 2.0, 3.0, 4.0}, tagged_alloc<double>{9});
    check(!(a == shorter), "mixed equality: prefix of equal values is not equality");
    check(!(a == longer), "mixed equality: extension of equal values is not equality");
  }

  // Scalars must match: a float vector and a double vector are not operands for each other.
  static_assert(
    !lam::linalg::compatible_vectors_c_weak<V<double, std::allocator<double>>, V<float, std::allocator<float>>>);
  static_assert(
    !requires(const V<double, std::allocator<double>>& x, const V<float, std::allocator<float>>& y) { x + y; });
  static_assert(!requires(const V<double, std::allocator<double>>& x, const V<float, std::allocator<float>>& y) {
    lam::linalg::dot(x, y);
  });

  // Vector-valued algorithms put the result on the left operand's allocator.
  {
    const std_v x({1.0, 0.0, 0.0});
    const tag_v y({0.0, 1.0, 0.0}, tagged_alloc<double>{3});

    auto cr = lam::linalg::cross(y, x);
    static_assert(std::same_as<decltype(cr), tag_v>);
    check(cr[2] == -1.0, "mixed cross values");
    check(cr.get_allocator().tag == 3, "mixed cross allocator");

    auto pr = lam::linalg::project(y, x);
    static_assert(std::same_as<decltype(pr), tag_v>);
    check(pr.get_allocator().tag == 3, "mixed project allocator");
    check(pr[0] == 0.0 && pr[1] == 0.0, "mixed project of orthogonal vectors is zero");

    auto rj = lam::linalg::reject(y, x);
    check(rj[1] == 1.0 && rj.get_allocator().tag == 3, "mixed reject");

    auto rf = lam::linalg::reflect(y, x);
    check(rf[1] == 1.0 && rf.get_allocator().tag == 3, "mixed reflect");

    auto lp = lam::linalg::lerp(y, x, 0.5);
    check(lp[0] == 0.5 && lp[1] == 0.5 && lp.get_allocator().tag == 3, "mixed lerp");

    auto mp = lam::linalg::midpoint(y, x);
    check(mp[0] == 0.5 && mp[1] == 0.5 && mp.get_allocator().tag == 3, "mixed midpoint");

    check(lam::linalg::is_orthogonal(x, y), "mixed is_orthogonal");
    check(!lam::linalg::is_parallel(x, y), "mixed is_parallel");

    const pmr_v z({0.0, 0.0, 1.0}, &pool);
    check(lam::linalg::triple_product(x, y, z) == 1.0, "mixed triple_product");
  }

  // The converting constructor is the only way to change allocator, and it always names one.
  {
    const tag_v moved(a, tagged_alloc<double>{42});
    check(moved.size() == 3 && moved[2] == 3.0, "converting ctor values");
    check(moved.get_allocator().tag == 42, "converting ctor allocator");
  }

  // Assignment across allocator types stays a compile error: it would have to pick one silently.
  static_assert(!std::assignable_from<std_v&, const tag_v&>);
  static_assert(!std::assignable_from<tag_v&, const pmr_v&>);

  std::println("{} mixed-allocator suite done", name);
}

// Mixing container types is the same rule, not a separate feature.
void cross_container()
{
  using v_t = lam::linalg::vector<double, std::allocator<double>>;
  using r_t = lam::linalg::r_vector<double, tagged_alloc<double>>;

  const v_t a({1.0, 2.0, 3.0});
  const r_t b({10.0, 20.0, 30.0}, tagged_alloc<double>{5});

  check(lam::linalg::dot(a, b) == 140.0, "cross-container dot");

  const v_t converted(b, std::allocator<double>{});
  check(converted[1] == 20.0, "cross-container converting ctor");

  auto s = lam::linalg::cross(a, b);
  static_assert(std::same_as<decltype(s), v_t>);
  check(s.size() == 3, "cross-container cross result type");
}

// The whole thing still works at compile time.
consteval bool constexpr_mixed()
{
  using std_v = lam::linalg::vector<double, std::allocator<double>>;
  using tag_v = lam::linalg::vector<double, tagged_alloc<double>>;

  std_v a({1.0, 2.0, 3.0});
  tag_v b({10.0, 20.0, 30.0}, tagged_alloc<double>{1});

  if (lam::linalg::dot(a, b) != 140.0)
    return false;

  auto s = a + b;
  if (s[0] != 11.0 || s[2] != 33.0)
    return false;

  tag_v converted(a, tagged_alloc<double>{2});
  if (converted[1] != 2.0 || converted.get_allocator().tag != 2)
    return false;

  a += b;
  return a[1] == 22.0;
}

static_assert(constexpr_mixed());

// Derived, not picked: a hand-picked 1e-10 sits below float's epsilon and breaks is_parallel.
template<typename T>
void angle_tolerance_suite(std::string_view name)
{
  using vec = lam::linalg::vector<T>;

  const vec a({T{1}, T{2}, T{3}});
  check(a.is_parallel(a), std::format("{}: a vector must be parallel to itself", name));
  check(a.is_parallel(vec({T{2}, T{4}, T{6}})), std::format("{}: exact scalar multiple is parallel", name));
  check(a.is_parallel(vec({T{-2}, T{-4}, T{-6}})), std::format("{}: negative multiple is parallel", name));
  check(!a.is_parallel(vec({T{1}, T{0}, T{0}})), std::format("{}: unrelated vectors are not parallel", name));

  const vec e1({T{1}, T{0}, T{0}});
  const vec e2({T{0}, T{1}, T{0}});
  check(e1.is_orthogonal(e2), std::format("{}: basis vectors are orthogonal", name));
  check(!e1.is_orthogonal(e1), std::format("{}: a vector is not orthogonal to itself", name));

  // The default must scale with the type, not be a constant.
  static_assert(lam::linalg::default_angle_tolerance<T> > std::numeric_limits<T>::epsilon());
  static_assert(lam::linalg::default_angle_tolerance<T> < T{1e-4});
}

int main()
{
  angle_tolerance_suite<float>("float");
  angle_tolerance_suite<double>("double");
  mixed_suite<lam::linalg::vector>("vector");
  mixed_suite<lam::linalg::r_vector>("r_vector");
  mixed_suite<lam::linalg::raw_vector>("raw_vector");
  cross_container();

  if (failures != 0)
  {
    std::println("{} failure(s)", failures);
    return 1;
  }
  std::println("\n✅ mixed-allocator tests passed");
  return 0;
}
