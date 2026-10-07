/*
 *  test_throwing_elements.cpp - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  test_throwing_elements is a c++ module
 */

import std;
import lam.linearalgebra;

using lam::linalg::storage_layout;

struct tracked_error : std::runtime_error
{
  tracked_error() : std::runtime_error("tracked") {}
};

// Scalar that counts live instances and throws once its construction budget runs out.
struct tracked
{
  static inline int live = 0;
  static inline int budget = std::numeric_limits<int>::max();

  double x = 0.0;

  static void spend()
  {
    if (budget == 0)
      throw tracked_error();
    --budget;
  }

  tracked(double d = 0.0) : x{d}
  {
    spend();
    ++live;
  }
  tracked(const tracked& o) : x{o.x}
  {
    spend();
    ++live;
  }
  tracked(tracked&& o) : x{o.x}
  {
    spend();
    ++live;
  }
  tracked& operator=(const tracked&) = default;
  tracked& operator=(tracked&&) = default;
  ~tracked() { --live; }

  tracked operator-() const { return tracked{-x}; }
  friend tracked operator+(const tracked& l, const tracked& r) { return tracked{l.x + r.x}; }
  friend tracked operator-(const tracked& l, const tracked& r) { return tracked{l.x - r.x}; }
  friend tracked operator*(const tracked& l, const tracked& r) { return tracked{l.x * r.x}; }
  friend bool operator==(const tracked& l, const tracked& r) { return l.x == r.x; }
};

// Counts live allocations.
struct counting_resource : std::pmr::memory_resource
{
  int live = 0;
  void* do_allocate(std::size_t n, std::size_t al) override
  {
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

using alloc_t = std::pmr::polymorphic_allocator<tracked>;

// Runs op with the throw placed at every construction in turn, until op completes.
// op(r1, r2, arm) builds its inputs, calls arm() to start the budget, then does the operation under test;
// it returns nonzero if a strong-guarantee check it performs fails.
int every_throw_point(auto op)
{
  for (int k = 0; k < 1000; ++k)
  {
    counting_resource r1;
    counting_resource r2;
    const int live_before = tracked::live;
    bool threw = false;
    int result = 0;
    try
    {
      result = op(r1, r2, [k] { tracked::budget = k; });
    }
    catch (const tracked_error&)
    {
      threw = true;
    }
    tracked::budget = std::numeric_limits<int>::max();
    if (result != 0)
      return result;
    if (r1.live != 0 || r2.live != 0)
      return 1000 + k;
    if (tracked::live != live_before)
      return 2000 + k;
    if (!threw)
      return k == 0 ? 3000 : 0; // k == 0 completing means the operation never constructed anything
  }
  return 4000;
}

// Assignment either completes or leaves the target exactly as it was.
template<typename V>
int unchanged(const V& target, std::initializer_list<double> expected)
{
  if (target.size() != expected.size())
    return 50;
  for (std::size_t i = 0; i < expected.size(); ++i)
    if (target[i].x != expected.begin()[i])
      return 51;
  return 0;
}

template<template<typename, typename> class V>
int vector_suite()
{
  using vec = V<tracked, alloc_t>;
  int r = 0;

  if ((r = every_throw_point([](auto& r1, auto&, auto arm) {
         arm();
         vec v(4, alloc_t{&r1});
         return 0;
       })))
    return 100 + r;

  if ((r = every_throw_point([](auto& r1, auto&, auto arm) {
         arm();
         vec v({tracked{1.0}, tracked{2.0}, tracked{3.0}}, alloc_t{&r1});
         return 0;
       })))
    return 200 + r;

  if ((r = every_throw_point([](auto& r1, auto&, auto arm) {
         std::vector<tracked> source{1.0, 2.0, 3.0};
         arm();
         vec v(source, alloc_t{&r1});
         return 0;
       })))
    return 300 + r;

  if ((r = every_throw_point([](auto& r1, auto&, auto arm) {
         vec source({1.0, 2.0, 3.0}, alloc_t{&r1});
         arm();
         vec copy(source);
         return 0;
       })))
    return 400 + r;

  if ((r = every_throw_point([](auto& r1, auto& r2, auto arm) {
         vec source({1.0, 2.0, 3.0}, alloc_t{&r1});
         arm();
         vec copy(source, alloc_t{&r2});
         return 0;
       })))
    return 500 + r;

  if ((r = every_throw_point([](auto& r1, auto& r2, auto arm) {
         vec target({7.0, 8.0}, alloc_t{&r2});
         vec source({1.0, 2.0, 3.0}, alloc_t{&r1});
         arm();
         try
         {
           target = source;
         }
         catch (const tracked_error&)
         {
           if (int u = unchanged(target, {7.0, 8.0}))
             return u;
           throw;
         }
         return unchanged(target, {1.0, 2.0, 3.0});
       })))
    return 600 + r;

  if ((r = every_throw_point([](auto& r1, auto& r2, auto arm) {
         vec target({7.0, 8.0}, alloc_t{&r2});
         vec source({1.0, 2.0, 3.0}, alloc_t{&r1});
         arm();
         try
         {
           target = std::move(source);
         }
         catch (const tracked_error&)
         {
           if (int u = unchanged(target, {7.0, 8.0}))
             return u;
           throw;
         }
         return unchanged(target, {1.0, 2.0, 3.0});
       })))
    return 700 + r;

  return 0;
}

int matrix_suite()
{
  using mat = lam::linalg::matrix<tracked, storage_layout::row_major, alloc_t>;
  int r = 0;

  if ((r = every_throw_point([](auto& r1, auto&, auto arm) {
         arm();
         mat m(2, 3, alloc_t{&r1});
         return 0;
       })))
    return 100 + r;

  if ((r = every_throw_point([](auto& r1, auto&, auto arm) {
         mat source(2, 3, alloc_t{&r1});
         arm();
         mat copy(source);
         return 0;
       })))
    return 200 + r;

  if ((r = every_throw_point([](auto& r1, auto& r2, auto arm) {
         mat target(2, 2, alloc_t{&r2});
         target[1, 1] = tracked{5.0};
         mat source(2, 3, alloc_t{&r1});
         arm();
         try
         {
           target = std::move(source);
         }
         catch (const tracked_error&)
         {
           if (target.rows() != 2 || target.cols() != 2 || target[1, 1].x != 5.0)
             return 52;
           throw;
         }
         return target.rows() == 2 && target.cols() == 3 ? 0 : 53;
       })))
    return 300 + r;

  if ((r = every_throw_point([](auto& r1, auto&, auto arm) {
         mat source(2, 3, alloc_t{&r1});
         arm();
         auto negated = -source;
         return 0;
       })))
    return 400 + r;

  return 0;
}

int failures = 0;

void check(int result, std::string_view what)
{
  if (result != 0)
  {
    std::println("FAILED {}: code {}", what, result);
    ++failures;
  }
  else
    std::println("{} passed", what);
}

int main()
{
  check(vector_suite<lam::linalg::vector>(), "vector, throwing element constructors");
  check(vector_suite<lam::linalg::r_vector>(), "r_vector, throwing element constructors");
  check(vector_suite<lam::linalg::raw_vector>(), "raw_vector, throwing element constructors");
  check(matrix_suite(), "matrix, throwing element constructors");

  if (failures != 0)
  {
    std::println("{} failure(s)", failures);
    return 1;
  }
  return 0;
}
