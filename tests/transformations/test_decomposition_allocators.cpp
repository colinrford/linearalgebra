/*
 *  test_decomposition_allocators.cpp - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  test_decomposition_allocators is a c++ module
 */

import std;
import lam.linearalgebra;

using namespace lam::linalg;

// Any global operator new during a decomposition means an allocation escaped the input's allocator.
std::atomic<long> heap_allocations{0};

// Every flavour of operator new is counted; only the aligned one needs the offset trick.
void* counted_allocate(std::size_t n)
{
  ++heap_allocations;
  if (void* p = std::malloc(n ? n : 1))
    return p;
  throw std::bad_alloc();
}

// Uncounted: the test's own resource allocates through this, and those are the allocations we want.
void* raw_aligned_alloc(std::size_t n, std::size_t align)
{
  align = std::max(align, alignof(void*));
  void* base = std::malloc(n + align + sizeof(void*));
  if (base == nullptr)
    throw std::bad_alloc();
  auto raw = reinterpret_cast<std::uintptr_t>(base) + sizeof(void*);
  auto aligned = (raw + align - 1) & ~(static_cast<std::uintptr_t>(align) - 1);
  std::memcpy(reinterpret_cast<void*>(aligned - sizeof(void*)), &base, sizeof(void*));
  return reinterpret_cast<void*>(aligned);
}

void* counted_allocate_aligned(std::size_t n, std::size_t align)
{
  ++heap_allocations;
  return raw_aligned_alloc(n, align);
}

void raw_aligned_free(void* p) noexcept
{
  if (p == nullptr)
    return;
  void* base = nullptr;
  std::memcpy(&base, static_cast<std::byte*>(p) - sizeof(void*), sizeof(void*));
  std::free(base);
}

void* operator new(std::size_t n) { return counted_allocate(n); }
void* operator new[](std::size_t n) { return counted_allocate(n); }
void* operator new(std::size_t n, std::align_val_t a) { return counted_allocate_aligned(n, std::to_underlying(a)); }
void* operator new[](std::size_t n, std::align_val_t a) { return counted_allocate_aligned(n, std::to_underlying(a)); }
void operator delete(void* p) noexcept { std::free(p); }
void operator delete[](void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }
void operator delete[](void* p, std::size_t) noexcept { std::free(p); }
void operator delete(void* p, std::align_val_t) noexcept { raw_aligned_free(p); }
void operator delete[](void* p, std::align_val_t) noexcept { raw_aligned_free(p); }
void operator delete(void* p, std::size_t, std::align_val_t) noexcept { raw_aligned_free(p); }
void operator delete[](void* p, std::size_t, std::align_val_t) noexcept { raw_aligned_free(p); }

// Counts live allocations.
struct counting_resource : std::pmr::memory_resource
{
  long live = 0;
  long total = 0;
  void* do_allocate(std::size_t n, std::size_t al) override
  {
    ++live;
    ++total;
    return raw_aligned_alloc(n, al);
  }
  void do_deallocate(void* p, std::size_t, std::size_t al) override
  {
    --live;
    raw_aligned_free(p);
  }
  bool do_is_equal(const std::pmr::memory_resource& o) const noexcept override { return this == &o; }
};

// Stateful non-pmr allocator; a default-constructed one throws, so a silent fallback fails loudly.
template<typename T>
struct arena_alloc
{
  using value_type = T;
  counting_resource* res = nullptr;

  arena_alloc() noexcept = default;
  explicit arena_alloc(counting_resource* r) noexcept : res{r} {}
  template<typename U>
  arena_alloc(const arena_alloc<U>& o) noexcept : res{o.res}
  {}

  T* allocate(std::size_t n)
  {
    if (res == nullptr)
      throw std::bad_alloc();
    return static_cast<T*>(res->allocate(n * sizeof(T), alignof(T)));
  }
  void deallocate(T* p, std::size_t n) { res->deallocate(p, n * sizeof(T), alignof(T)); }

  template<typename U>
  friend bool operator==(const arena_alloc& a, const arena_alloc<U>& b) noexcept
  { return a.res == b.res; }
};

int failures = 0;

void check(bool ok, std::string_view what)
{
  if (!ok)
  {
    std::println("FAILED {}", what);
    ++failures;
  }
}

// Runs f with heap allocations counted; returns how many happened.
long heap_during(auto&& f)
{
  long before = heap_allocations.load();
  f();
  return heap_allocations.load() - before;
}

template<typename Alloc>
bool on(const Alloc& a, counting_resource* r)
{
  if constexpr (requires { a.resource(); })
    return a.resource() == r;
  else
    return a.res == r;
}

template<typename Alloc, typename T = double>
using mat_t = matrix<T, storage_layout::row_major, typename std::allocator_traits<Alloc>::template rebind_alloc<T>>;

template<template<typename> class AllocFor>
void run_suite(std::string_view name, auto make_alloc)
{
  counting_resource r;
  {
    using A = AllocFor<double>;
    using C = AllocFor<std::complex<double>>;
    A alloc = make_alloc(&r);
    C calloc_ = C(alloc);

    mat_t<A> a(3, 3, alloc);
    a[0, 0] = 4.0;
    a[0, 1] = 1.0;
    a[1, 0] = 1.0;
    a[1, 1] = 3.0;
    a[2, 2] = 2.0;
    mat_t<A> b(3, 3, alloc);
    b[0, 0] = b[1, 1] = b[2, 2] = 1.0;
    vector<double, A> rhs({1.0, 2.0, 3.0}, alloc);

    mat_t<C, std::complex<double>> ca(2, 2, calloc_);
    ca[0, 0] = {2.0, 0.0};
    ca[0, 1] = {0.0, 1.0};
    ca[1, 0] = {0.0, -1.0};
    ca[1, 1] = {3.0, 0.0};
    mat_t<C, std::complex<double>> cb(2, 2, calloc_);
    cb[0, 0] = cb[1, 1] = {1.0, 0.0};

    std::optional<svd_result<double, A>> s;
    check(heap_during([&] { s = svd(a); }) == 0, std::format("{}: svd allocated on the heap", name));
    check(s.has_value() && s->s.size() == 3, std::format("{}: svd failed", name));
    if (s)
    {
      check(on(s->u.get_allocator(), &r) && on(s->vt.get_allocator(), &r) && on(s->s.get_allocator(), &r),
            std::format("{}: svd result not on the input's resource", name));
      check(std::abs(s->s[0] - (7.0 + std::sqrt(5.0)) / 2.0) < 1e-12 &&
              std::abs(s->s[1] - (7.0 - std::sqrt(5.0)) / 2.0) < 1e-12 && std::abs(s->s[2] - 2.0) < 1e-12,
            std::format("{}: svd singular values", name));
      auto diag = svd_check(a, *s);
      check(diag.reconstruction_error < 1e-12, std::format("{}: svd reconstruction", name));
    }

    std::optional<mat_t<A>> p;
    check(heap_during([&] { p = pinv(a); }) == 0, std::format("{}: pinv allocated on the heap", name));
    check(on(p->get_allocator(), &r), std::format("{}: pinv result not on the input's resource", name));
    auto identity = a * *p;
    check(std::abs(identity[0, 0] - 1.0) < 1e-12 && std::abs(identity[0, 1]) < 1e-12,
          std::format("{}: pinv is not an inverse", name));

    std::optional<vector<double, A>> x;
    check(heap_during([&] { x = lstsq(a, rhs); }) == 0, std::format("{}: lstsq allocated on the heap", name));
    check(on(x->get_allocator(), &r), std::format("{}: lstsq result not on the input's resource", name));

    std::size_t rank = 0;
    check(heap_during([&] { rank = matrix_rank(a); }) == 0, std::format("{}: matrix_rank allocated on the heap", name));
    check(rank == 3, std::format("{}: matrix_rank", name));

    std::optional<eigen_result<double, A>> e;
    check(heap_during([&] { e = symmetric_eigen(a); }) == 0,
          std::format("{}: symmetric_eigen allocated on the heap", name));
    check(e && e->success && std::abs(e->eigenvalues[0] - 2.0) < 1e-12,
          std::format("{}: symmetric_eigen values", name));
    check(e && on(e->eigenvalues.get_allocator(), &r) && on(e->eigenvectors.get_allocator(), &r),
          std::format("{}: symmetric_eigen result not on the input's resource", name));

    std::optional<gep_result<double, A>> g;
    check(heap_during([&] { g = solve_gep(a, b); }) == 0, std::format("{}: solve_gep allocated on the heap", name));
    check(g && g->success && std::abs(g->eigenvalues[0] - 2.0) < 1e-12, std::format("{}: solve_gep values", name));
    check(g && on(g->eigenvalues.get_allocator(), &r) && on(g->eigenvectors.get_allocator(), &r),
          std::format("{}: solve_gep result not on the input's resource", name));

    // Past n = 100 solve_gep switches to the divide-and-conquer LAPACK path, with its own workspaces.
    constexpr std::size_t big = 120;
    mat_t<A> big_a(big, big, alloc);
    mat_t<A> big_b(big, big, alloc);
    for (std::size_t i = 0; i < big; ++i)
    {
      big_a[i, i] = static_cast<double>(i + 1);
      big_b[i, i] = 1.0;
    }
    std::optional<gep_result<double, A>> bg;
    check(heap_during([&] { bg = solve_gep(big_a, big_b); }) == 0,
          std::format("{}: divide-and-conquer solve_gep allocated on the heap", name));
    check(bg && bg->success && std::abs(bg->eigenvalues[0] - 1.0) < 1e-12 && on(bg->eigenvalues.get_allocator(), &r),
          std::format("{}: divide-and-conquer solve_gep result", name));

    std::optional<gep_result<std::complex<double>, C>> cg;
    check(heap_during([&] { cg = solve_gep(ca, cb); }) == 0,
          std::format("{}: complex solve_gep allocated on the heap", name));
    check(cg && cg->success && on(cg->eigenvectors.get_allocator(), &r),
          std::format("{}: complex solve_gep result not on the input's resource", name));

    std::optional<eigen_result<std::complex<double>, C>> ce;
    check(heap_during([&] { ce = general_complex_eigen(ca); }) == 0,
          std::format("{}: general_complex_eigen allocated on the heap", name));
    check(ce && ce->success && on(ce->eigenvectors.get_allocator(), &r),
          std::format("{}: general_complex_eigen result not on the input's resource", name));

    // float and complex SVD take different LAPACK branches, each with their own workspaces.
    using F = AllocFor<float>;
    F falloc = F(alloc);
    mat_t<F, float> fa(3, 3, falloc);
    fa[0, 0] = 4.0F;
    fa[1, 1] = 3.0F;
    fa[2, 2] = 2.0F;
    std::optional<svd_result<float, F>> fs;
    check(heap_during([&] { fs = svd(fa); }) == 0, std::format("{}: float svd allocated on the heap", name));
    check(fs && on(fs->u.get_allocator(), &r) && on(fs->s.get_allocator(), &r),
          std::format("{}: float svd result not on the input's resource", name));
    check(fs && std::abs(fs->s[0] - 4.0F) < 1e-5F, std::format("{}: float svd singular values", name));

    using CF = AllocFor<std::complex<float>>;
    CF cfalloc = CF(alloc);
    mat_t<CF, std::complex<float>> cfa(2, 2, cfalloc);
    cfa[0, 0] = {4.0F, 0.0F};
    cfa[1, 1] = {2.0F, 0.0F};
    std::optional<svd_result<std::complex<float>, CF>> cfs;
    check(heap_during([&] { cfs = svd(cfa); }) == 0, std::format("{}: complex<float> svd allocated on the heap", name));
    check(cfs && on(cfs->u.get_allocator(), &r) && on(cfs->s.get_allocator(), &r),
          std::format("{}: complex<float> svd result not on the input's resource", name));
    check(cfs && std::abs(cfs->s[0] - 4.0F) < 1e-5F, std::format("{}: complex<float> svd singular values", name));

    std::optional<svd_result<std::complex<double>, C>> cs;
    check(heap_during([&] { cs = svd(ca); }) == 0, std::format("{}: complex svd allocated on the heap", name));
    check(cs && on(cs->u.get_allocator(), &r) && on(cs->s.get_allocator(), &r),
          std::format("{}: complex svd result not on the input's resource", name));

    // Failure paths return empty results that still carry the input's allocator.
    mat_t<A> rect(2, 3, alloc);
    auto bad = symmetric_eigen(rect);
    check(!bad.success && bad.eigenvalues.size() == 0 && on(bad.eigenvectors.get_allocator(), &r),
          std::format("{}: symmetric_eigen failure result", name));
    auto bad_gep = solve_gep(rect, rect);
    check(!bad_gep.success && on(bad_gep.eigenvalues.get_allocator(), &r),
          std::format("{}: solve_gep failure result", name));
    auto bad_complex = general_complex_eigen(mat_t<C, std::complex<double>>(2, 3, calloc_));
    check(!bad_complex.success && on(bad_complex.eigenvalues.get_allocator(), &r),
          std::format("{}: general_complex_eigen failure result", name));
  }
  check(r.live == 0, std::format("{}: leaked from the resource", name));
  check(r.total > 0, std::format("{}: nothing was allocated on the resource", name));
  std::println("{} done", name);
}

template<typename T>
using pmr_alloc = std::pmr::polymorphic_allocator<T>;

int main()
{
  std::pmr::set_default_resource(std::pmr::null_memory_resource());

  // Default allocators still produce default-allocator results.
  static_assert(std::same_as<decltype(svd(matrix<double>(2, 2)))::value_type, svd_result<double>>);
  static_assert(std::same_as<decltype(symmetric_eigen(matrix<double>(2, 2))), eigen_result<double>>);
  static_assert(std::same_as<decltype(pinv(matrix<double>(2, 2))), matrix<double>>);

  run_suite<pmr_alloc>("polymorphic_allocator", [](counting_resource* r) { return pmr_alloc<double>{r}; });
  run_suite<arena_alloc>("arena_alloc", [](counting_resource* r) { return arena_alloc<double>{r}; });

  if (failures != 0)
  {
    std::println("{} failure(s)", failures);
    return 1;
  }
  std::println("All decomposition allocator tests passed");
  return 0;
}
