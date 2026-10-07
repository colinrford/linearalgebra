/*
 *  test_allocator_propagation.cpp - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  test_allocator_propagation is a c++ module
 */

import std;
import lam.linearalgebra;
import lam.concepts;

// Counts live allocations; throws bad_alloc once `allow` runs out.
struct counting_resource : std::pmr::memory_resource
{
  std::size_t live = 0;
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

// Scalar whose constructor throws once its budget runs out.
struct fussy
{
  static inline int budget = 0;
  double x = 0.0;
  fussy(double d = 0.0) : x{d}
  {
    if (budget-- == 0)
      throw std::runtime_error("fussy");
  }
  fussy(const fussy& o) : fussy(o.x) {}
  fussy& operator=(const fussy&) = default;
  fussy operator-() const { return fussy{-x}; }
  friend fussy operator+(const fussy& l, const fussy& r) { return fussy{l.x + r.x}; }
  friend fussy operator-(const fussy& l, const fussy& r) { return fussy{l.x - r.x}; }
  friend fussy operator*(const fussy& l, const fussy& r) { return fussy{l.x * r.x}; }
  friend bool operator==(const fussy&, const fussy&) = default;
};

// Helper to check conditions
void check(bool condition, std::string_view msg)
{
  if (!condition)
  {
    throw std::runtime_error(std::string(msg));
  }
}

int main()
{
  using namespace lam::linalg;

  // Create a buffer for PMR
  std::array<std::byte, 4096> buffer;
  std::pmr::monotonic_buffer_resource pool{buffer.data(), buffer.size()};
  std::pmr::polymorphic_allocator<double> alloc{&pool};

  using vec_t = vector<double, std::pmr::polymorphic_allocator<double>>;

  // Test 1: Construction with allocator
  {
    vec_t v({1.0, 2.0, 3.0}, alloc);
    check(v.get_allocator() == alloc, "Construction: Allocator not stored correctly");
    std::println("Construction passed");
  }

  // Test 2: Copy propagation
  {
    vec_t v1({1.0, 2.0, 3.0}, alloc);
    vec_t v2 = v1; // Copy constructor
    check(v2.get_allocator() == alloc, "Copy: Allocator not propagated");
    std::println("Copy propagation passed");
  }

  // Test 3: Allocator-extended copy
  {
    vec_t v1({1.0, 2.0, 3.0}, alloc);
    // Create a different allocator (null resource for testing difference)
    std::pmr::polymorphic_allocator<double> alloc2{};
    vec_t v2(v1, alloc2);
    check(v2.get_allocator() == alloc2, "Extended copy: Allocator not used");
    check(v2.get_allocator() != v1.get_allocator(), "Extended copy: Allocator should be different");
    std::println("Allocator-extended copy passed");
  }

  // Test 4: Operator+ propagation
  {
    vec_t v1({1.0, 2.0, 3.0}, alloc);
    vec_t v2({4.0, 5.0, 6.0}, alloc);
    vec_t v3 = v1 + v2;
    check(v3.get_allocator() == alloc, "Operator+: Allocator not propagated");
    std::println("Operator+ passed");
  }

  // Test 5: Operator* propagation
  {
    vec_t v1({1.0, 2.0, 3.0}, alloc);
    vec_t v2 = v1 * 2.0;
    check(v2.get_allocator() == alloc, "Operator*: Allocator not propagated");
    vec_t v3 = 2.0 * v1;
    check(v3.get_allocator() == alloc, "Operator* (scalar first): Allocator not propagated");
    std::println("Operator* passed");
  }

  // Test 6: Member functions (cross, unit, etc.)
  {
    vec_t v1({1.0, 0.0, 0.0}, alloc);
    vec_t v2({0.0, 1.0, 0.0}, alloc);

    vec_t v_cross = v1.cross(v2);
    check(v_cross.get_allocator() == alloc, "cross(): Allocator not propagated");

    vec_t v_unit = v1.unit();
    check(v_unit.get_allocator() == alloc, "unit(): Allocator not propagated");

    vec_t v_project = v1.project(v2);
    check(v_project.get_allocator() == alloc, "project(): Allocator not propagated"); // Should take from *this (v1)

    vec_t v_lerp = v1.lerp(v2, 0.5);
    check(v_lerp.get_allocator() == alloc, "lerp(): Allocator not propagated");

    std::println("Member functions passed");
  }

  // Test 7: Propagation vs Non-Propagation (Comparison)
  {
    vec_t v1({10.0, 20.0}, alloc);

    // Case A: Propagation (Default Copy Constructor)
    vec_t v_prop = v1;
    check(v_prop.get_allocator() == alloc, "Case A: Default copy should propagate allocator");
    check(v_prop.get_allocator() == v1.get_allocator(), "Case A: Allocators should match");

    // Case B: Non-Propagation (Explicitly using default allocator)
    std::pmr::polymorphic_allocator<double> default_alloc{}; // Uses default resource (new/delete)
    vec_t v_noprop(v1, default_alloc);

    check(v_noprop.get_allocator() == default_alloc, "Case B: Explicit copy should use provided allocator");
    check(v_noprop.get_allocator() != alloc, "Case B: Should NOT use source allocator");
    check(v_noprop.get_allocator().resource() != alloc.resource(),
          "Case B: Resources should differ (Default vs Monotonic)");

    std::println("Propagation vs Non-Propagation verified");
  }

  // Test 8: Matrix results adopt the operand's resource, not the default one
  {
    using mat_t = matrix<double, storage_layout::row_major, std::pmr::polymorphic_allocator<double>>;
    mat_t a(2, 3, alloc);
    mat_t b(3, 2, alloc);
    vec_t v({1.0, 2.0, 3.0}, alloc);

    check((a * b).get_allocator() == alloc, "matrix*matrix: Allocator not propagated");
    check(a.transpose().get_allocator() == alloc, "transpose(): Allocator not propagated");
    check((a * v).get_allocator() == alloc, "matrix*vector: Allocator not propagated");

    std::println("Matrix result propagation passed");
  }

  // Test 9: raw_vector stays consistent when allocation or construction throws
  {
    using raw_t = raw_vector<double, std::pmr::polymorphic_allocator<double>>;
    counting_resource once;
    once.allow = 1;
    std::pmr::polymorphic_allocator<double> throwing{&once};

    // Unequal resources + different size forces a fresh allocation on the target.
    raw_t target({1.0, 2.0}, throwing);
    raw_t source({1.0, 2.0, 3.0}, alloc);
    bool threw = false;
    try
    {
      target = std::move(source);
    }
    catch (const std::bad_alloc&)
    {
      threw = true;
    }
    check(threw, "raw_vector move-assign: expected bad_alloc from exhausted resource");
    check(target.size() == 2 && target[0] == 1.0 && target[1] == 2.0, "raw_vector move-assign: strong guarantee lost");

    std::println("raw_vector throw safety passed");
  }

  // Test 10: raw_vector frees partially-built storage when an element ctor throws
  {
    counting_resource res;
    std::pmr::polymorphic_allocator<fussy> falloc{&res};
    fussy::budget = 2; // third element construction throws
    bool threw = false;
    try
    {
      raw_vector<fussy, std::pmr::polymorphic_allocator<fussy>> fv(4, falloc);
    }
    catch (const std::runtime_error&)
    {
      threw = true;
    }
    check(threw, "raw_vector ctor: expected element ctor to throw");
    check(res.live == 0, "raw_vector ctor: storage leaked after element ctor threw");

    std::println("raw_vector ctor rollback passed");
  }

  // Test 11: LU results, and unary minus, stay on the input's resource
  {
    using mat_t = matrix<double, storage_layout::row_major, std::pmr::polymorphic_allocator<double>>;
    mat_t a(2, 2, alloc);
    a[0, 0] = 4.0;
    a[0, 1] = 3.0;
    a[1, 0] = 6.0;
    a[1, 1] = 3.0;
    vec_t b({1.0, 2.0}, alloc);
    vector<double> b_default({1.0, 2.0});

    check((-a).get_allocator() == alloc, "matrix unary minus: Allocator not propagated");
    check(crout_lu(a)->lu_decomp.get_allocator() == alloc, "crout_lu: Allocator not propagated");
    check(crout_lu(a)->permutations.get_allocator() == alloc, "crout_lu perms: Allocator not propagated");
    check(doolittle_lu(a)->lu_decomp.get_allocator() == alloc, "doolittle_lu: Allocator not propagated");
    check(crout_lu_solve(a, b).get_allocator() == alloc, "crout_lu_solve: Allocator not propagated");
    check(crout_lu_solve(a, b_default).get_allocator() == alloc,
          "crout_lu_solve (default b): Allocator not propagated");
    check(crout_lu_inv(a).get_allocator() == alloc, "crout_lu_inv: Allocator not propagated");

    std::println("LU propagation passed");
  }

  // Test 12: assignment works with pmr, keeping each target's own resource
  {
    using mat_t = matrix<double, storage_layout::row_major, std::pmr::polymorphic_allocator<double>>;
    std::array<std::byte, 4096> other_buffer;
    std::pmr::monotonic_buffer_resource other_pool{other_buffer.data(), other_buffer.size()};
    std::pmr::polymorphic_allocator<double> other_alloc{&other_pool};

    vec_t a({1.0, 2.0, 3.0}, alloc);
    vec_t b({4.0, 5.0, 6.0}, other_alloc);

    b = a; // copy-assign across unequal resources
    check(b[0] == 1.0 && b[2] == 3.0, "pmr copy-assign: values not copied");
    check(b.get_allocator() == other_alloc, "pmr copy-assign: target must keep its own allocator");

    vec_t c({7.0, 8.0}, other_alloc);
    c = vec_t({9.0, 9.0, 9.0}, alloc); // move-assign, different size, unequal resources
    check(c.size() == 3 && c[0] == 9.0, "pmr move-assign: values not moved");
    check(c.get_allocator() == other_alloc, "pmr move-assign: target must keep its own allocator");

    mat_t m(2, 2, alloc);
    mat_t n(2, 2, other_alloc);
    m[0, 0] = 2.0;
    n = m;
    check(n[0, 0] == 2.0, "pmr matrix copy-assign: values not copied");
    check(n.get_allocator() == other_alloc, "pmr matrix assign: target must keep its own allocator");

    mat_t p(2, 2, alloc);
    p[0, 1] = 5.0;
    n = p * m;
    check(n.get_allocator() == other_alloc, "pmr matrix assign from expression: allocator changed");

    std::println("pmr assignment passed");
  }

  // Test 13: the pmr aliases name the same types, and deduce their resource from the allocator
  {
    namespace pmr = lam::linalg::pmr;

    static_assert(std::same_as<pmr::vector<double>, vector<double, std::pmr::polymorphic_allocator<double>>>);
    static_assert(
      std::same_as<pmr::r_vector<double>, lam::linalg::r_vector<double, std::pmr::polymorphic_allocator<double>>>);
    static_assert(
      std::same_as<pmr::raw_vector<double>, lam::linalg::raw_vector<double, std::pmr::polymorphic_allocator<double>>>);
    static_assert(std::same_as<pmr::matrix<double>,
                               matrix<double, storage_layout::row_major, std::pmr::polymorphic_allocator<double>>>);
    static_assert(std::same_as<pmr::matrix<double, storage_layout::col_major>,
                               matrix<double, storage_layout::col_major, std::pmr::polymorphic_allocator<double>>>);

    // A resource converts to the allocator, so no allocator type is ever spelled at the call site.
    std::array<std::byte, 4096> alias_buffer;
    std::pmr::monotonic_buffer_resource alias_pool{alias_buffer.data(), alias_buffer.size()};

    pmr::vector<double> v({1.0, 2.0, 3.0}, &alias_pool);
    check(v.size() == 3 && v[2] == 3.0, "pmr::vector alias: values wrong");
    check(v.get_allocator().resource() == &alias_pool, "pmr::vector alias: resource not adopted");

    pmr::matrix<double> m(2, 2, &alias_pool);
    m[0, 0] = 4.0;
    check(m.get_allocator().resource() == &alias_pool, "pmr::matrix alias: resource not adopted");

    // Results of operations stay on the same resource.
    auto sum = m + m;
    check(sum[0, 0] == 8.0, "pmr::matrix alias: sum wrong");
    check(sum.get_allocator().resource() == &alias_pool, "pmr::matrix alias: result left the resource");

    std::println("pmr aliases passed");
  }

  std::println("\n✅ All allocator propagation tests passed!");
  return 0;
}
