/*
 *  vectorspace_algorithms.cppm - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  vectorspace_algorithms is a c++ module
 */

export module lam.linearalgebra:vectorspace.algorithms;

import std;
import lam.concepts;

import :vectorspace.exceptions;
import :config;

namespace lam::linalg
{

#ifdef LAM_USE_BLAS
extern "C"
{
  double cblas_ddot(const int N, const double* X, const int incX, const double* Y, const int incY);
  float cblas_sdot(const int N, const float* X, const int incX, const float* Y, const int incY);
  double cblas_dnrm2(const int N, const double* X, const int incX);
  float cblas_snrm2(const int N, const float* X, const int incX);
}

template<typename T>
struct blas_vector_dispatcher;

template<>
struct blas_vector_dispatcher<double>
{
  static double dot(int N, const double* X, int incX, const double* Y, int incY)
  {
    return cblas_ddot(N, X, incX, Y, incY);
  }
  static double nrm2(int N, const double* X, int incX) { return cblas_dnrm2(N, X, incX); }
};

template<>
struct blas_vector_dispatcher<float>
{
  static float dot(int N, const float* X, int incX, const float* Y, int incY)
  {
    return cblas_sdot(N, X, incX, Y, incY);
  }
  static float nrm2(int N, const float* X, int incX) { return cblas_snrm2(N, X, incX); }
};
#endif

constexpr auto sqrt_helper(auto x)
{
  using T = decltype(x);
  if (x < T{0})
    throw vector_exception::non_pos(); // Handle domain error
  if (x == T{0})
    return T{0};
  T g = x;
  for (int i = 0; i < 20; ++i)
  { // 20 iterations is plenty for double precision
    g = (g + x / g) / 2;
  }
  return g;
}

export template<lam::concepts::experimental::vector_c_weak V>
constexpr typename V::scalar_type dot(const V& a, const V& b)
{
  using T = typename V::scalar_type;
  typename V::size_type min_dim = std::min(a.size(), b.size());

  auto generic_impl = [&]() {
    T sum = T{0};
    for (typename V::size_type i = 0; i < min_dim; ++i)
    {
      sum += a[i] * b[i];
    }
    return sum;
  };

  if consteval
  {
    return generic_impl();
  }
  else
  {
    if constexpr (config::use_blas && (std::is_same_v<T, double> || std::is_same_v<T, float>))
    {
#ifdef LAM_USE_BLAS
      if constexpr (std::ranges::contiguous_range<V>)
      {
        if (min_dim == 0)
          return T{0};
        const T* pa = std::to_address(std::ranges::begin(a));
        const T* pb = std::to_address(std::ranges::begin(b));
        // Note: Assuming stride 1 for contiguous vectors
        return blas_vector_dispatcher<T>::dot(static_cast<int>(min_dim), pa, 1, pb, 1);
      }
      else
      {
        return generic_impl();
      }
#else
      return generic_impl();
#endif
    }
    else
    {
      return generic_impl();
    }
  }
}

export template<lam::concepts::experimental::vector_c_weak V>
constexpr typename V::scalar_type norm2(const V& v)
{
  using T = typename V::scalar_type;

  if consteval
    return dot(v, v);
  else
  {
    if constexpr (config::use_blas && (std::is_same_v<T, double> || std::is_same_v<T, float>))
    {
#ifdef LAM_USE_BLAS
      if constexpr (std::ranges::contiguous_range<V>)
      {
        if (v.size() == 0)
          return T{0};
        const T* pv = std::to_address(std::ranges::begin(v));
        // Use nrm2 squared? No cblas_dnrm2 returns the norm (sqrt(dot(v,v))).
        // So for norm2 (squared norm), we might just stick to dot(v,v) or square the result of nrm2.
        // However, dot(v,v) via BLAS (cblas_ddot) is efficient.
        // Using cblas_dnrm2 involves a square root internally, and then we square it? Wasteful.
        // Better to use cblas_ddot(v, v) for norm2.
        // So actually, just calling dot(v, v) is fine because dot is already optimized!
        return dot(v, v);
      }
      else
        return dot(v, v);
#else
      return dot(v, v);
#endif
    }
    else
      return dot(v, v);
  }
}

export 
template<lam::concepts::experimental::vector_c_weak V>
constexpr typename V::scalar_type norm(const V& v)
{
  using T = typename V::scalar_type;
  if consteval
    return sqrt_helper(norm2(v));
  else
  {
    if constexpr (config::use_blas && (std::is_same_v<T, double> || std::is_same_v<T, float>))
    {
#ifdef LAM_USE_BLAS
      if constexpr (std::ranges::contiguous_range<V>)
      {
        if (v.size() == 0)
          return T{0};
        const T* pv = std::to_address(std::ranges::begin(v));
        return blas_vector_dispatcher<T>::nrm2(static_cast<int>(v.size()), pv, 1);
      }
#endif
    }
    return std::sqrt(norm2(v));
  }
}

export 
template<std::ranges::range R1, std::ranges::range R2>
  requires requires(R1 r1, R2 r2, std::size_t i) {
    r1[i];
    r2[i];
    { r1.size() } -> std::convertible_to<std::size_t>;
  }
constexpr auto dot_range(const R1& r1, const R2& r2)
{
  using T = std::common_type_t<std::ranges::range_value_t<R1>, std::ranges::range_value_t<R2>>;

  if consteval
  { // Use iota-based approach for constexpr (zip has limitations 1/8/26)
    auto n = std::min(r1.size(), r2.size());
    auto indices = std::views::iota(std::size_t{0}, n);
    return std::ranges::fold_left(indices, T{0}, [&](T sum, auto i) { return sum + r1[i] * r2[i]; });
  }
  else
  { // Use zip-based approach for runtime performance
    return std::ranges::fold_left(std::views::zip(r1, r2), T{0},
                                  [](T sum, const auto& pair) { return sum + std::get<0>(pair) * std::get<1>(pair); });
  }
}

export 
template<std::ranges::range R>
constexpr auto norm2_range(const R& r)
{ return dot_range(r, r); }

export 
template<std::ranges::range R>
constexpr auto norm_range(const R& r)
{
  if consteval
    return sqrt_helper(norm2_range(r));
  else
    return std::sqrt(norm2_range(r));
}

export 
template<lam::concepts::experimental::vector_c_weak V>
constexpr typename V::scalar_type distance(const V& a, const V& b)
{
  using T = typename V::scalar_type;
  typename V::size_type min_dim = std::min(a.size(), b.size());
  T sum = T{0};

  for (typename V::size_type i = 0; i < min_dim; ++i)
  {
    T d = static_cast<T>(a[i]) - static_cast<T>(b[i]);
    sum += d * d;
  }

  if consteval
    return sqrt_helper(sum);
  else
    return std::sqrt(sum);
}

export 
template<std::ranges::range R1, std::ranges::range R2>
  requires requires(R1 r1, R2 r2, std::size_t i) {
    r1[i];
    r2[i];
    { r1.size() } -> std::convertible_to<std::size_t>;
  }
constexpr auto distance_range(const R1& r1, const R2& r2)
{
  using T = std::common_type_t<std::ranges::range_value_t<R1>, std::ranges::range_value_t<R2>>;

  if consteval
  { // Use iota-based approach for constexpr (zip has limitations)
    auto n = std::min(r1.size(), r2.size());
    auto indices = std::views::iota(std::size_t{0}, n);
    auto sum_sq = std::ranges::fold_left(indices, T{0}, [&](T sum, auto i) {
      auto d = r1[i] - r2[i];
      return sum + d * d;
    });
    return sqrt_helper(sum_sq);
  }
  else
  { // Use zip-based approach for runtime performance
    auto sum_sq = std::ranges::fold_left(std::views::zip(r1, r2), T{0}, [](T sum, const auto& pair) {
      auto d = std::get<0>(pair) - std::get<1>(pair);
      return sum + d * d;
    });
    return std::sqrt(sum_sq);
  }
}

export 
template<lam::concepts::experimental::vector_c_weak V>
constexpr typename V::scalar_type angle(const V& a, const V& b)
{
  using T = typename V::scalar_type;
  T n1 = norm(a);
  T n2 = norm(b);
  if (n1 == T{0} || n2 == T{0})
    throw vector_exception::div_by_zero();
  T cos_theta = dot(a, b) / (n1 * n2);
  if (cos_theta > T{1})
    cos_theta = T{1};
  if (cos_theta < T{-1})
    cos_theta = T{-1};

  if consteval
  {
    T x = cos_theta;
    T x2 = x * x;
    T x3 = x2 * x;
    T x5 = x3 * x2;
    T x7 = x5 * x2;
    T asin_approx = x + x3 / T{6} + T{3} * x5 / T{40} + T{15} * x7 / T{336};
    constexpr T pi_over_2 = T{1.5707963267948966};
    return pi_over_2 - asin_approx;
  }
  else
  {
    return std::acos(cos_theta);
  }
}

export 
template<lam::concepts::experimental::vector_c_weak V>
constexpr V cross(const V& a, const V& b)
{
  if (a.size() != 3 || b.size() != 3)
    throw vector_exception::cross_undef();

  if constexpr (requires { a.get_allocator(); })
  {
    V res(3, a.get_allocator());
    res[0] = a[1] * b[2] - a[2] * b[1];
    res[1] = a[2] * b[0] - a[0] * b[2];
    res[2] = a[0] * b[1] - a[1] * b[0];
    return res;
  }
  else
  {
    V res(3);
    res[0] = a[1] * b[2] - a[2] * b[1];
    res[1] = a[2] * b[0] - a[0] * b[2];
    res[2] = a[0] * b[1] - a[1] * b[0];
    return res;
  }
}

export 
template<lam::concepts::experimental::vector_c_weak V>
constexpr V lerp(const V& a, const V& b, typename V::scalar_type t)
{
  using T = typename V::scalar_type;
  return (T{1} - t) * a + t * b;
}

export 
template<lam::concepts::experimental::vector_c_weak V>
constexpr V project(const V& v, const V& onto)
{
  using T = typename V::scalar_type;
  T on2 = norm2(onto);
  if (on2 == T{0})
    throw vector_exception::div_by_zero();
  T s = dot(v, onto) / on2;
  return s * onto;
}

export 
template<lam::concepts::experimental::vector_c_weak V>
constexpr V reject(const V& v, const V& from)
{ return v - project(v, from); }

export 
template<lam::concepts::experimental::vector_c_weak V>
constexpr V reflect(const V& v, const V& normal)
{
  // r = v - 2(v·n̂)n̂  where n̂ is the unit normal
  using T = typename V::scalar_type;
  T n = norm(normal);
  if (n == T{0})
    throw vector_exception::div_by_zero();

  // Scale normal to unit length
  V n_unit = normal;
  n_unit = (T{1} / n) * n_unit;

  T scale = T{2} * dot(v, n_unit);
  return v - scale * n_unit;
}

export 
template<lam::concepts::experimental::vector_c_weak V>
constexpr V unit(const V& v)
{
  using T = typename V::scalar_type;
  T n = norm(v);
  if (n == T{0})
    throw vector_exception::div_by_zero();
  return (T{1} / n) * v;
}

export 
template<lam::concepts::experimental::vector_c_weak V>
constexpr V midpoint(const V& a, const V& b)
{
  using T = typename V::scalar_type;
  return (T{0.5}) * (a + b);
}

export 
template<lam::concepts::experimental::vector_c_weak V>
constexpr bool is_parallel(const V& a, const V& b, typename V::scalar_type tolerance = typename V::scalar_type{1e-10})
{
  using T = typename V::scalar_type;
  T n1 = norm(a);
  T n2 = norm(b);
  if (n1 == T{0} || n2 == T{0})
    return true;

  T cos_theta = dot(a, b) / (n1 * n2);
  T abs_cos = cos_theta < T{0} ? -cos_theta : cos_theta;
  return (T{1} - abs_cos) < tolerance;
}

export 
template<lam::concepts::experimental::vector_c_weak V>
constexpr bool is_orthogonal(const V& a, const V& b, typename V::scalar_type tolerance = typename V::scalar_type{1e-10})
{
  using T = typename V::scalar_type;
  T n1 = norm(a);
  T n2 = norm(b);
  if (n1 == T{0} || n2 == T{0})
    return true;

  T cos_theta = dot(a, b) / (n1 * n2);
  T abs_cos = cos_theta < T{0} ? -cos_theta : cos_theta;
  return abs_cos < tolerance;
}

export 
template<lam::concepts::experimental::vector_c_weak V>
constexpr typename V::scalar_type triple_product(const V& a, const V& b, const V& c)
{ return dot(a, cross(b, c)); }

} // namespace lam::linalg
