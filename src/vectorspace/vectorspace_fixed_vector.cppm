/*
 *  vectorspace_fixed_vector.cppm - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  fixed_vector: A stack-allocated, fixed-size vector type.
 *  Satisfies vector_c_weak concept, enabling use with lam::linalg algorithms.
 *  Designed for constexpr contexts and graphics/physics applications.
 *
 *  NOTE: This partition is NOT automatically exported.
 *  To enable, add to vectorspace.cppm:
 *    export import :vectorspace.fixed_vector;
 */

module;

export module lam.linearalgebra:vectorspace.fixed_vector;

import std;
import lam.concepts;

namespace lam::linalg
{

// ============================================================================
// fixed_vector<T, N> - Stack-allocated, fixed-size vector
// ============================================================================
export 
template<typename T, std::size_t N>
  requires lam::concepts::experimental::ring_element_c_weak<T>
struct fixed_vector
{
  using scalar_type = T;
  using size_type = std::size_t;
  using value_type = T;

  std::array<T, N> data{};
  // ========== Constructors ==========
  constexpr fixed_vector() = default;

  // Zero constructor: fixed_vector(0) creates zero vector (required for vector_c_weak)
  constexpr explicit fixed_vector(int zero_val) noexcept
  {
    T zero_elem = static_cast<T>(zero_val);
    for (size_type i = 0; i < N; ++i)
      data[i] = zero_elem;
  }

  // Static zero factory (required for vector_c_weak)
  static constexpr fixed_vector zero() noexcept
  {
    return fixed_vector(0);
  }

  // Variadic constructor: fixed_vector<double, 3>{1.0, 2.0, 3.0}
  template<typename... Args>
    requires (sizeof...(Args) == N) && (std::convertible_to<Args, T> && ...)
  constexpr fixed_vector(Args... args) noexcept
    : data{static_cast<T>(args)...}
  {}
  // From std::array
  constexpr explicit fixed_vector(const std::array<T, N>& arr) noexcept
    : data{arr}
  {}
  // ========== Size (required for vector_c_weak) ==========
  static constexpr size_type static_size = N;
  constexpr size_type size() const noexcept { return N; }
  // ========== Element Access ==========
  constexpr scalar_type& operator[](size_type i) noexcept { return data[i]; }
  constexpr const scalar_type& operator[](size_type i) const noexcept { return data[i]; }

  constexpr scalar_type& at(size_type i)
  {
    if (i >= N) throw std::out_of_range("fixed_vector index out of bounds");
    return data[i];
  }
  constexpr const scalar_type& at(size_type i) const
  {
    if (i >= N) throw std::out_of_range("fixed_vector index out of bounds");
    return data[i];
  }
  // ========== Named Accessors (conditional on size) ==========
  constexpr scalar_type& x() noexcept requires (N >= 1) { return data[0]; }
  constexpr scalar_type& y() noexcept requires (N >= 2) { return data[1]; }
  constexpr scalar_type& z() noexcept requires (N >= 3) { return data[2]; }
  constexpr scalar_type& w() noexcept requires (N >= 4) { return data[3]; }

  constexpr const scalar_type& x() const noexcept requires (N >= 1) { return data[0]; }
  constexpr const scalar_type& y() const noexcept requires (N >= 2) { return data[1]; }
  constexpr const scalar_type& z() const noexcept requires (N >= 3) { return data[2]; }
  constexpr const scalar_type& w() const noexcept requires (N >= 4) { return data[3]; }
  // ========== Iterators ==========
  constexpr scalar_type* begin() noexcept { return data.data(); }
  constexpr scalar_type* end() noexcept { return data.data() + N; }
  constexpr const scalar_type* begin() const noexcept { return data.data(); }
  constexpr const scalar_type* end() const noexcept { return data.data() + N; }
  constexpr const scalar_type* cbegin() const noexcept { return data.data(); }
  constexpr const scalar_type* cend() const noexcept { return data.data() + N; }
  // ========== Arithmetic (required for module_element_c_weak) ==========
  constexpr fixed_vector operator-() const noexcept
  {
    fixed_vector result;
    for (size_type i = 0; i < N; ++i)
      result.data[i] = -data[i];
    return result;
  }

  constexpr fixed_vector& operator+=(const fixed_vector& other) noexcept
  {
    for (size_type i = 0; i < N; ++i)
      data[i] += other.data[i];
    return *this;
  }

  constexpr fixed_vector& operator-=(const fixed_vector& other) noexcept
  {
    for (size_type i = 0; i < N; ++i)
      data[i] -= other.data[i];
    return *this;
  }

  constexpr fixed_vector& operator*=(const scalar_type& scalar) noexcept
  {
    for (size_type i = 0; i < N; ++i)
      data[i] *= scalar;
    return *this;
  }

  constexpr fixed_vector& operator/=(const scalar_type& scalar) noexcept
  {
    for (size_type i = 0; i < N; ++i)
      data[i] /= scalar;
    return *this;
  }
  // ========== Comparison ==========
  constexpr bool operator==(const fixed_vector&) const noexcept = default;
  // ========== Raw Data Access (for GPU upload) ==========
  constexpr const scalar_type* raw() const noexcept { return data.data(); }
  constexpr scalar_type* raw() noexcept { return data.data(); }
};

// ============================================================================
// Free Function Operators
// ============================================================================
export 
template<typename T, std::size_t N>
constexpr fixed_vector<T, N> operator+(const fixed_vector<T, N>& a,
                                        const fixed_vector<T, N>& b) noexcept
{
  fixed_vector<T, N> result = a;
  result += b;
  return result;
}

export 
template<typename T, std::size_t N>
constexpr fixed_vector<T, N> operator-(const fixed_vector<T, N>& a,
                                        const fixed_vector<T, N>& b) noexcept
{
  fixed_vector<T, N> result = a;
  result -= b;
  return result;
}

export 
template<typename T, std::size_t N>
constexpr fixed_vector<T, N> operator*(const T& scalar,
                                        const fixed_vector<T, N>& v) noexcept
{
  fixed_vector<T, N> result = v;
  result *= scalar;
  return result;
}

export 
template<typename T, std::size_t N>
constexpr fixed_vector<T, N> operator*(const fixed_vector<T, N>& v,
                                        const T& scalar) noexcept
{
  return scalar * v;
}

export 
template<typename T, std::size_t N>
constexpr fixed_vector<T, N> operator/(const fixed_vector<T, N>& v,
                                        const T& scalar) noexcept
{
  fixed_vector<T, N> result = v;
  result /= scalar;
  return result;
}

// ============================================================================
// Concept Verification
// ============================================================================

static_assert(lam::concepts::experimental::vector_c_weak<fixed_vector<double, 3>>,
              "fixed_vector<double, 3> must satisfy vector_c_weak");
static_assert(lam::concepts::experimental::vector_c_weak<fixed_vector<float, 4>>,
              "fixed_vector<float, 4> must satisfy vector_c_weak");
static_assert(lam::concepts::experimental::vector_c_weak<fixed_vector<int, 2>>,
              "fixed_vector<int, 2> must satisfy vector_c_weak");

// ============================================================================
// Common Type Aliases
// ============================================================================

export using vec2f = fixed_vector<float, 2>;
export using vec3f = fixed_vector<float, 3>;
export using vec4f = fixed_vector<float, 4>;

export using vec2d = fixed_vector<double, 2>;
export using vec3d = fixed_vector<double, 3>;
export using vec4d = fixed_vector<double, 4>;

export using vec2i = fixed_vector<int, 2>;
export using vec3i = fixed_vector<int, 3>;
export using vec4i = fixed_vector<int, 4>;

} // namespace lam::linalg

// ============================================================================
// std::formatter specialization
// ============================================================================
export 
template<typename T, std::size_t N>
struct std::formatter<lam::linalg::fixed_vector<T, N>>
{
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  auto format(const lam::linalg::fixed_vector<T, N>& v,
              std::format_context& ctx) const
  {
    auto out = ctx.out();
    std::format_to(out, "<");
    for (std::size_t i = 0; i < N; ++i)
    {
      if (i > 0) 
        std::format_to(out, ", ");
      std::format_to(out, "{}", v[i]);
    }
    return std::format_to(out, ">");
  }
};
