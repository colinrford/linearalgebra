/*
 *  vectorspace_r_vector.cppm - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  vectorspace_r_vector is a c++ module
 */

export module lam.linearalgebra:vectorspace.r_vector;

import std;
import lam.concepts;

import :vectorspace.exceptions;
import :vectorspace.algorithms;

namespace lam::linalg
{

export template<typename T, typename Alloc = std::allocator<T>>
  requires lam::concepts::experimental::ring_element_c_weak<T>
class r_vector
{
public:
  using scalar_type = T;
  using allocator_type = Alloc;
  using size_type = std::size_t;
  using pointer = scalar_type*;
  using const_pointer = const scalar_type*;

private:
  // A throwing ~T() would std::terminate during cleanup; reject it at compile time.
  static_assert(std::is_nothrow_destructible_v<T>,
                "lam::linalg::r_vector requires T to have a non-throwing destructor");

  struct Deleter
  {
    Alloc alloc;
    size_type n;
    constexpr Deleter(const Alloc& a = Alloc(), size_type s = 0) : alloc(a), n(s) {}
    constexpr void operator()(pointer p)
    {
      if (p)
      {
        for (size_type i = 0; i < n; ++i)
        {
          std::allocator_traits<Alloc>::destroy(alloc, p + i);
        }
        std::allocator_traits<Alloc>::deallocate(alloc, p, n);
      }
    }
  };

  [[no_unique_address]] Alloc m_alloc;
  size_type m_dimension;
  std::unique_ptr<scalar_type[], Deleter> arrow;

public:
  constexpr r_vector() : m_alloc{}, m_dimension{0}, arrow{nullptr, Deleter{}} {}

  constexpr explicit r_vector(const allocator_type& alloc) noexcept
    : m_alloc{alloc}, m_dimension{0}, arrow{nullptr, Deleter{alloc, 0}}
  {}

  constexpr explicit r_vector(size_type dim)
    requires std::default_initializable<allocator_type>
    : r_vector(dim, allocator_type())
  {}

  constexpr r_vector(size_type dim, const allocator_type& alloc)
    : m_alloc{alloc}, m_dimension{dim}, arrow{nullptr, Deleter{alloc, dim}}
  {
    if (dim == 0)
      throw vector_exception::non_pos();
    arrow.reset(detail::build_storage(
      m_alloc, dim, [&](pointer p, size_type) { std::allocator_traits<Alloc>::construct(m_alloc, p); }));
  }

  constexpr r_vector(std::initializer_list<scalar_type> il)
    requires std::default_initializable<allocator_type>
    : r_vector(il, allocator_type())
  {}

  constexpr r_vector(std::initializer_list<scalar_type> il, const allocator_type& alloc)
    : m_alloc{alloc}, m_dimension{il.size()}, arrow{nullptr, Deleter{alloc, il.size()}}
  {
    if (il.size() == 0)
      throw vector_exception::non_pos();
    arrow.reset(detail::build_storage(m_alloc, m_dimension, [&](pointer p, size_type i) {
      std::allocator_traits<Alloc>::construct(m_alloc, p, il.begin()[i]);
    }));
  }

  constexpr r_vector(const std::vector<scalar_type>& v)
    requires std::default_initializable<allocator_type>
    : r_vector(v, allocator_type())
  {}

  constexpr r_vector(const std::vector<scalar_type>& v, const allocator_type& alloc)
    : m_alloc{alloc}, m_dimension{v.size()}, arrow{nullptr, Deleter{alloc, v.size()}}
  {
    if (v.size() == 0)
      throw vector_exception::non_pos();
    arrow.reset(detail::build_storage(m_alloc, m_dimension, [&](pointer p, size_type i) {
      std::allocator_traits<Alloc>::construct(m_alloc, p, v[i]);
    }));
  }

  constexpr r_vector(r_vector&& other) noexcept
    : m_alloc{other.m_alloc}, m_dimension{std::exchange(other.m_dimension, 0)}, arrow{std::move(other.arrow)}
  { other.arrow.get_deleter().n = 0; }

  // Steal storage only if the allocator propagates or compares equal; otherwise move into our own storage.
  constexpr r_vector& operator=(r_vector&& other) noexcept(detail::steals_on_move_assignment<Alloc>)
  {
    if (this == &other)
      return *this;

    if (detail::propagates_on_move_assignment<Alloc> || m_alloc == other.m_alloc)
    {
      arrow.reset();
      if constexpr (detail::propagates_on_move_assignment<Alloc>)
      {
        m_alloc = std::move(other.m_alloc);
        arrow.get_deleter().alloc = m_alloc;
      }
      arrow.get_deleter().n = other.m_dimension;
      arrow.reset(other.arrow.release());
      m_dimension = other.m_dimension;
      other.arrow.get_deleter().n = 0;
      other.m_dimension = 0;
    }
    else if (other.m_dimension == 0)
    {
      arrow.reset();
      arrow.get_deleter().n = 0;
      m_dimension = 0;
    }
    else if (m_dimension == other.m_dimension)
    {
      for (size_type i = 0; i < m_dimension; ++i)
        arrow[i] = std::move(other.arrow[i]);
    }
    else
    {
      pointer fresh = detail::build_storage(m_alloc, other.m_dimension, [&](pointer p, size_type i) {
        std::allocator_traits<Alloc>::construct(m_alloc, p, std::move(other.arrow[i]));
      });
      arrow.reset();
      arrow.get_deleter().n = other.m_dimension;
      arrow.reset(fresh);
      m_dimension = other.m_dimension;
    }
    return *this;
  }

  constexpr r_vector(const r_vector& other)
    : m_alloc{other.m_alloc}, m_dimension{other.m_dimension}, arrow{nullptr, Deleter{m_alloc, other.m_dimension}}
  {
    if (m_dimension > 0)
      arrow.reset(detail::build_storage(m_alloc, m_dimension, [&](pointer p, size_type i) {
        std::allocator_traits<Alloc>::construct(m_alloc, p, other.arrow[i]);
      }));
  }

  constexpr r_vector(const r_vector& other, const allocator_type& alloc)
    : m_alloc{alloc}, m_dimension{other.m_dimension}, arrow{nullptr, Deleter{alloc, other.m_dimension}}
  {
    if (m_dimension > 0)
      arrow.reset(detail::build_storage(m_alloc, m_dimension, [&](pointer p, size_type i) {
        std::allocator_traits<Alloc>::construct(m_alloc, p, other.arrow[i]);
      }));
  }

  // Any vector over the same scalars; the allocator is named here, never inherited from the source.
  template<typename Other>
    requires(!std::same_as<Other, r_vector>) && compatible_vectors_c_weak<r_vector, Other>
  constexpr r_vector(const Other& other, const allocator_type& alloc)
    : m_alloc{alloc}, m_dimension{other.size()}, arrow{nullptr, Deleter{alloc, other.size()}}
  {
    if (m_dimension > 0)
    {
      auto src = other.as_span();
      arrow.reset(detail::build_storage(m_alloc, m_dimension, [&](pointer p, size_type i) {
        std::allocator_traits<Alloc>::construct(m_alloc, p, src[i]);
      }));
    }
  }

  constexpr r_vector& operator=(const r_vector& other)
  {
    if (this != &other)
    {
      r_vector tmp(other, m_alloc);
      *this = std::move(tmp);
    }
    return *this;
  }

  constexpr scalar_type& operator[](size_type index)
  {
    if (index >= m_dimension)
      throw vector_exception::out_of_bounds();
    return arrow[index];
  }

  constexpr const scalar_type& operator[](size_type index) const
  {
    if (index >= m_dimension)
      throw vector_exception::out_of_bounds();
    return arrow[index];
  }

  constexpr size_type size() const noexcept { return m_dimension; }
  constexpr size_type dim() const noexcept { return m_dimension; }
  constexpr size_type dimension() const noexcept { return m_dimension; }
  constexpr size_type length() const noexcept { return m_dimension; }

  constexpr allocator_type get_allocator() const noexcept { return m_alloc; }

  // Throws when the allocators differ and do not propagate; swap never copies elements.
  constexpr void swap(r_vector& other) noexcept(detail::steals_on_swap<Alloc>)
  {
    if constexpr (detail::propagates_on_swap<Alloc>)
    {
      using std::swap;
      swap(m_alloc, other.m_alloc);
    }
    else if (!(m_alloc == other.m_alloc))
      throw vector_exception::unequal_allocators();
    std::swap(m_dimension, other.m_dimension);
    arrow.swap(other.arrow);
  }
  friend constexpr void swap(r_vector& a, r_vector& b) noexcept(noexcept(a.swap(b))) { a.swap(b); }

  // Degenerate additive identity. A zero with a dimension needs one: zero_like(v).
  static constexpr r_vector zero()
    requires std::default_initializable<allocator_type>
  { return r_vector{}; }

  // Generic Delegates
  template<typename Other>
    requires compatible_vectors_c_weak<r_vector, Other>
  constexpr scalar_type dot(const Other& other) const
  { return lam::linalg::dot(*this, other); }
  constexpr scalar_type norm2() const { return lam::linalg::norm2(*this); }
  constexpr scalar_type norm() const { return lam::linalg::norm(*this); }

  constexpr scalar_type pnorm(scalar_type p) const
  {
    auto abs_pow = std::views::transform(as_span(), [p](auto val) { return std::pow(std::abs(val), p); });
    scalar_type sum = std::ranges::fold_left(abs_pow, scalar_type{0}, std::plus<scalar_type>{});
    return std::pow(sum, scalar_type{1} / p);
  }

  template<typename Other>
    requires compatible_vectors_c_weak<r_vector, Other>
  constexpr r_vector cross(const Other& other) const
  { return lam::linalg::cross(*this, other); }
  constexpr r_vector unit() const { return lam::linalg::unit(*this); }

  template<typename Other>
    requires compatible_vectors_c_weak<r_vector, Other>
  constexpr bool operator==(const Other& other) const
  {
    // ranges::equal already reports unequal lengths as not-equal.
    return std::ranges::equal(as_span(), other.as_span());
  }

  constexpr r_vector operator-() const
  {
    r_vector res(*this, m_alloc);
    std::ranges::transform(res.as_span(), res.begin(), [](auto val) { return -val; });
    return res;
  }

  template<typename Other>
    requires compatible_vectors_c_weak<r_vector, Other>
  constexpr r_vector& operator+=(const Other& other)
  {
    if (m_dimension != other.size())
      throw vector_exception::plus_equals_unequal_dim();

    for (auto&& [a, b] : std::views::zip(as_span(), other.as_span()))
    {
      a += b;
    }
    return *this;
  }

  constexpr pointer begin() noexcept { return arrow.get(); }
  constexpr pointer end() noexcept { return arrow.get() + m_dimension; }
  constexpr const_pointer begin() const noexcept { return arrow.get(); }
  constexpr const_pointer end() const noexcept { return arrow.get() + m_dimension; }
  constexpr const_pointer cbegin() const noexcept { return arrow.get(); }
  constexpr const_pointer cend() const noexcept { return arrow.get() + m_dimension; }

  constexpr std::span<scalar_type> as_span() noexcept { return {arrow.get(), m_dimension}; }
  constexpr std::span<const scalar_type> as_span() const noexcept { return {arrow.get(), m_dimension}; }

  template<typename Other>
    requires compatible_vectors_c_weak<r_vector, Other>
  constexpr r_vector& operator-=(const Other& other)
  {
    if (m_dimension != other.size())
      throw vector_exception::plus_equals_unequal_dim();

    for (auto&& [a, b] : std::views::zip(as_span(), other.as_span()))
    {
      a -= b;
    }
    return *this;
  }

  constexpr r_vector& operator*=(const scalar_type& s)
  {
    std::ranges::for_each(as_span(), [s](auto& val) { val *= s; });
    return *this;
  }

  constexpr r_vector& operator/=(const scalar_type& s)
  {
    if (s == scalar_type{0})
      throw vector_exception::div_by_zero();
    scalar_type inv_s = scalar_type{1} / s;
    std::ranges::for_each(as_span(), [inv_s](auto& val) { val *= inv_s; });
    return *this;
  }

  constexpr r_vector& normalize()
  {
    scalar_type n = norm();
    if (n == scalar_type{0})
      throw vector_exception::div_by_zero();
    scalar_type inv_n = scalar_type{1} / n;
    operator*=(inv_n);
    return *this;
  }

  // Generic Delegates
  template<typename Other>
    requires compatible_vectors_c_weak<r_vector, Other>
  constexpr scalar_type angle(const Other& other) const
  { return lam::linalg::angle(*this, other); }
  template<typename Other>
    requires compatible_vectors_c_weak<r_vector, Other>
  constexpr r_vector project(const Other& onto) const
  { return lam::linalg::project(*this, onto); }
  template<typename Other>
    requires compatible_vectors_c_weak<r_vector, Other>
  constexpr r_vector reject(const Other& from) const
  { return lam::linalg::reject(*this, from); }
  template<typename Other>
    requires compatible_vectors_c_weak<r_vector, Other>
  constexpr r_vector reflect(const Other& normal) const
  { return lam::linalg::reflect(*this, normal); }
  template<typename Other>
    requires compatible_vectors_c_weak<r_vector, Other>
  constexpr r_vector lerp(const Other& other, scalar_type t) const
  { return lam::linalg::lerp(*this, other, t); }
  template<typename Other>
    requires compatible_vectors_c_weak<r_vector, Other>
  constexpr scalar_type distance(const Other& other) const
  { return lam::linalg::distance(*this, other); }
  template<typename Other>
    requires compatible_vectors_c_weak<r_vector, Other>
  constexpr bool is_parallel(const Other& other, scalar_type tolerance = default_angle_tolerance<scalar_type>) const
  { return lam::linalg::is_parallel(*this, other, tolerance); }
  template<typename Other>
    requires compatible_vectors_c_weak<r_vector, Other>
  constexpr bool is_orthogonal(const Other& other, scalar_type tolerance = default_angle_tolerance<scalar_type>) const
  { return lam::linalg::is_orthogonal(*this, other, tolerance); }
  template<typename B, typename C>
    requires compatible_vectors_c_weak<r_vector, B> && compatible_vectors_c_weak<r_vector, C>
  constexpr scalar_type triple_product(const B& b, const C& c) const
  { return lam::linalg::triple_product(*this, b, c); }
};

// Operators
// Operands may carry different allocators; the result takes the LEFT one's.
export template<typename T, typename A1, typename A2>
constexpr r_vector<T, A1> operator+(const r_vector<T, A1>& a, const r_vector<T, A2>& b)
{
  r_vector<T, A1> res(a, a.get_allocator());
  res += b;
  return res;
}

// Operands may carry different allocators; the result takes the LEFT one's.
export template<typename T, typename A1, typename A2>
constexpr r_vector<T, A1> operator-(const r_vector<T, A1>& a, const r_vector<T, A2>& b)
{
  r_vector<T, A1> res(a, a.get_allocator());
  res -= b;
  return res;
}

export template<typename T, typename Alloc>
constexpr r_vector<T, Alloc> operator*(const T& s, const r_vector<T, Alloc>& v)
{
  r_vector<T, Alloc> res(v, v.get_allocator());
  res *= s;
  return res;
}

export template<typename T, typename Alloc>
constexpr r_vector<T, Alloc> operator*(const r_vector<T, Alloc>& v, const T& s)
{ return s * v; }

export template<typename T, typename Alloc>
constexpr r_vector<T, Alloc> operator/(const r_vector<T, Alloc>& v, const T& s)
{
  r_vector<T, Alloc> res(v, v.get_allocator());
  res /= s;
  return res;
}

// Spares the caller from naming an allocator type: a memory_resource* converts to the allocator.
export namespace pmr
{
template<typename T>
using r_vector = lam::linalg::r_vector<T, std::pmr::polymorphic_allocator<T>>;
} // namespace pmr

} // namespace lam::linalg
