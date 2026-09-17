/*
 *  bench_support.cppm - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  Module linalg.bench.support: shared scaffolding for the google/benchmark
 *  micro suite. Imports std + lam.linearalgebra ONLY -- deliberately NO
 *  google/benchmark dependency here. The DoNotOptimize/ClobberMemory barrier
 *  lives in the .cpp benchmark TUs, where <benchmark/benchmark.h> is already
 *  included; this module carries only dependency-clean, reusable pieces.
 */

export module linalg.bench.support;

import std;
import lam.linearalgebra;

export namespace linalg::bench
{

// Application allocation profiles the micro benchmarks model.
enum class profile
{
  games,
  graphics3d,
  ml
};

constexpr std::string_view name_of(profile p) noexcept
{
  switch (p)
  {
    case profile::games:
      return "games";
    case profile::graphics3d:
      return "graphics3d";
    case profile::ml:
      return "ml";
  }
  return "unknown";
}

// Deterministic (fixed-seed) value source so runs are reproducible and no
// benchmark accidentally measures RNG. splitmix64 -> double in [0, 1).
class deterministic_source
{
  std::uint64_t state_;

public:
  constexpr explicit deterministic_source(std::uint64_t seed = 0x9E3779B97F4A7C15ull) noexcept : state_{seed} {}

  constexpr double next() noexcept
  {
    state_ += 0x9E3779B97F4A7C15ull;
    std::uint64_t z = state_;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    z ^= (z >> 31);
    return static_cast<double>(z >> 11) * (1.0 / 9007199254740992.0);
  }
};

// Fill anything with a 1-arg operator[] (e.g. lam::linalg::vector) of known length.
// Constrained on 1-arg subscript so it can't collide with the matrix overload
// below (both are 3-arg-callable once the seed defaults).
template<typename Indexable>
  requires requires(Indexable& x, std::size_t i) { x[i]; }
constexpr void fill(Indexable& x, std::size_t n, std::uint64_t seed = 0x1234ull)
{
  deterministic_source src{seed};
  for (std::size_t i = 0; i < n; ++i)
    x[i] = src.next();
}

// Fill anything with a 2-arg operator[] (e.g. lam::linalg::matrix) of known shape.
template<typename Matrix>
  requires requires(Matrix& m, std::size_t i, std::size_t j) { m[i, j]; }
constexpr void fill(Matrix& m, std::size_t rows, std::size_t cols, std::uint64_t seed = 0x1234ull)
{
  deterministic_source src{seed};
  for (std::size_t i = 0; i < rows; ++i)
    for (std::size_t j = 0; j < cols; ++j)
      m[i, j] = src.next();
}

// RAII monotonic arena for the games profile: owns a fixed byte buffer + resource,
// hands out a polymorphic_allocator, and reset() rewinds it (the per-frame move).
template<typename T>
class monotonic_arena
{
  std::vector<std::byte> buffer_; // declared before mr_: initialized first
  std::pmr::monotonic_buffer_resource mr_;

public:
  explicit monotonic_arena(std::size_t bytes) : buffer_(bytes), mr_{buffer_.data(), buffer_.size()} {}

  std::pmr::polymorphic_allocator<T> allocator() noexcept { return {&mr_}; }

  // Rewind to the owned buffer's start -- models the per-frame arena reset.
  void reset() noexcept { mr_.release(); }
};

} // namespace linalg::bench
