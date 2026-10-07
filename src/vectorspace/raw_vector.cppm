/*
 *  vectorspace_raw_vector.cppm - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  vectorspace_raw_vector is a c++ module
 */

export module lam.linearalgebra:vectorspace.raw_vector;

import std;
import lam.concepts;

import :vectorspace.exceptions;
import :vectorspace.algorithms;

namespace lam::linalg
{

export template<typename T, typename Alloc = std::allocator<T>>
  requires lam::concepts::experimental::ring_element_c_weak<T>
class raw_vector
{
public:
  using scalar_type = T;
  using allocator_type = Alloc;
  using size_type = std::size_t;
  using pointer = scalar_type*;
  using const_pointer = const scalar_type*;

  // Throwing destructors would std::terminate during cleanup; reject at compile time.
  static_assert(std::is_nothrow_destructible_v<T>,
                "lam::linalg::raw_vector requires T to have a non-throwing destructor");

private:
  [[no_unique_address]] Alloc m_alloc;
  size_type m_dimension;
  pointer m_data; // owning; nullptr iff m_dimension == 0

  // Destroy (non-trivial T only) then deallocate our storage. Idempotent.
  constexpr void deallocate_storage() noexcept
  {
    if (m_data)
    {
      if constexpr (!std::is_trivially_destructible_v<scalar_type>)
        for (size_type i = 0; i < m_dimension; ++i)
          std::allocator_traits<Alloc>::destroy(m_alloc, m_data + i);
      std::allocator_traits<Alloc>::deallocate(m_alloc, m_data, m_dimension);
      m_data = nullptr;
    }
  }

public:
  constexpr raw_vector() noexcept : m_alloc{}, m_dimension{0}, m_data{nullptr} {}

  constexpr explicit raw_vector(const allocator_type& alloc) noexcept : m_alloc{alloc}, m_dimension{0}, m_data{nullptr}
  {}

  constexpr explicit raw_vector(size_type dim)
    requires std::default_initializable<allocator_type>
    : raw_vector(dim, allocator_type())
  {}

  constexpr raw_vector(size_type dim, const allocator_type& alloc) : m_alloc{alloc}, m_dimension{dim}, m_data{nullptr}
  {
    if (dim == 0)
      throw vector_exception::non_pos();
    m_data = detail::build_storage(m_alloc, dim,
                                   [&](pointer p, size_type) { std::allocator_traits<Alloc>::construct(m_alloc, p); });
  }

  constexpr raw_vector(std::initializer_list<scalar_type> il)
    requires std::default_initializable<allocator_type>
    : raw_vector(il, allocator_type())
  {}

  constexpr raw_vector(std::initializer_list<scalar_type> il, const allocator_type& alloc)
    : m_alloc{alloc}, m_dimension{il.size()}, m_data{nullptr}
  {
    if (il.size() == 0)
      throw vector_exception::non_pos();
    m_data = detail::build_storage(m_alloc, m_dimension, [&](pointer p, size_type i) {
      std::allocator_traits<Alloc>::construct(m_alloc, p, il.begin()[i]);
    });
  }

  constexpr raw_vector(const std::vector<scalar_type>& v)
    requires std::default_initializable<allocator_type>
    : raw_vector(v, allocator_type())
  {}

  constexpr raw_vector(const std::vector<scalar_type>& v, const allocator_type& alloc)
    : m_alloc{alloc}, m_dimension{v.size()}, m_data{nullptr}
  {
    if (v.size() == 0)
      throw vector_exception::non_pos();
    m_data = detail::build_storage(
      m_alloc, m_dimension, [&](pointer p, size_type i) { std::allocator_traits<Alloc>::construct(m_alloc, p, v[i]); });
  }

  constexpr ~raw_vector() { deallocate_storage(); }

  // Steal the pointer and null the source; a raw owning pointer needs this by hand.
  constexpr raw_vector(raw_vector&& other) noexcept
    : m_alloc{std::move(other.m_alloc)}, m_dimension{other.m_dimension}, m_data{other.m_data}
  {
    other.m_data = nullptr;
    other.m_dimension = 0;
  }

  // Steal only under POCMA or an equal allocator; otherwise move element-wise into our own storage.
  constexpr raw_vector& operator=(raw_vector&& other) noexcept(detail::steals_on_move_assignment<Alloc>)
  {
    if (this == &other)
      return *this;

    if (detail::propagates_on_move_assignment<Alloc> || m_alloc == other.m_alloc)
    {
      deallocate_storage();
      if constexpr (detail::propagates_on_move_assignment<Alloc>)
        m_alloc = std::move(other.m_alloc);
      m_data = other.m_data;
      m_dimension = other.m_dimension;
      other.m_data = nullptr;
      other.m_dimension = 0;
    }
    else if (m_dimension == other.m_dimension)
    {
      for (size_type i = 0; i < m_dimension; ++i)
        m_data[i] = std::move(other.m_data[i]);
    }
    else
    {
      // Build first, then swap in: a throw leaves *this untouched.
      pointer fresh =
        other.m_dimension
          ? detail::build_storage(m_alloc, other.m_dimension,
                                  [&](pointer p, size_type i) {
                                    std::allocator_traits<Alloc>::construct(m_alloc, p, std::move(other.m_data[i]));
                                  })
          : nullptr;
      deallocate_storage();
      m_data = fresh;
      m_dimension = other.m_dimension;
    }
    return *this;
  }

  // Propagates the source allocator: a deliberate deviation from std::pmr::vector.
  constexpr raw_vector(const raw_vector& other)
    : m_alloc{other.m_alloc}, m_dimension{other.m_dimension}, m_data{nullptr}
  {
    if (m_dimension > 0)
      m_data = detail::build_storage(m_alloc, m_dimension, [&](pointer p, size_type i) {
        std::allocator_traits<Alloc>::construct(m_alloc, p, other[i]);
      });
  }

  constexpr raw_vector(const raw_vector& other, const allocator_type& alloc)
    : m_alloc{alloc}, m_dimension{other.m_dimension}, m_data{nullptr}
  {
    if (m_dimension > 0)
      m_data = detail::build_storage(m_alloc, m_dimension, [&](pointer p, size_type i) {
        std::allocator_traits<Alloc>::construct(m_alloc, p, other[i]);
      });
  }

  // Any vector over the same scalars; the allocator is named here, never inherited from the source.
  template<typename Other>
    requires(!std::same_as<Other, raw_vector>) && compatible_vectors_c_weak<raw_vector, Other>
  constexpr raw_vector(const Other& other, const allocator_type& alloc)
    : m_alloc{alloc}, m_dimension{other.size()}, m_data{nullptr}
  {
    if (m_dimension > 0)
    {
      auto src = other.as_span();
      m_data = detail::build_storage(m_alloc, m_dimension, [&](pointer p, size_type i) {
        std::allocator_traits<Alloc>::construct(m_alloc, p, src[i]);
      });
    }
  }

  constexpr raw_vector& operator=(const raw_vector& other)
  {
    if (this != &other)
    {
      raw_vector tmp(other, m_alloc);
      *this = std::move(tmp);
    }
    return *this;
  }

  constexpr scalar_type& operator[](size_type index)
  {
    if (index >= m_dimension)
      throw vector_exception::out_of_bounds();
    return m_data[index];
  }

  constexpr const scalar_type& operator[](size_type index) const
  {
    if (index >= m_dimension)
      throw vector_exception::out_of_bounds();
    return m_data[index];
  }

  constexpr size_type size() const noexcept { return m_dimension; }
  constexpr size_type dim() const noexcept { return m_dimension; }
  constexpr size_type dimension() const noexcept { return m_dimension; }
  constexpr size_type length() const noexcept { return m_dimension; }

  constexpr allocator_type get_allocator() const noexcept { return m_alloc; }

  // Throws when the allocators differ and do not propagate; swap never copies elements.
  constexpr void swap(raw_vector& other) noexcept(detail::steals_on_swap<Alloc>)
  {
    if constexpr (detail::propagates_on_swap<Alloc>)
    {
      using std::swap;
      swap(m_alloc, other.m_alloc);
    }
    else if (!(m_alloc == other.m_alloc))
      throw vector_exception::unequal_allocators();
    std::swap(m_dimension, other.m_dimension);
    std::swap(m_data, other.m_data);
  }
  friend constexpr void swap(raw_vector& a, raw_vector& b) noexcept(noexcept(a.swap(b))) { a.swap(b); }

  // Degenerate additive identity. A zero with a dimension needs one: zero_like(v).
  static constexpr raw_vector zero()
    requires std::default_initializable<allocator_type>
  { return raw_vector{}; }

  template<typename Other>
    requires compatible_vectors_c_weak<raw_vector, Other>
  constexpr scalar_type dot(const Other& other) const
  { return lam::linalg::dot(*this, other); }
  constexpr scalar_type norm2() const { return lam::linalg::norm2(*this); }
  constexpr scalar_type norm() const { return lam::linalg::norm(*this); }

  constexpr scalar_type pnorm(scalar_type p) const
  {
    scalar_type sum = scalar_type{0};
    for (const auto& val : *this)
    {
      sum += std::pow(std::abs(val), p);
    }
    return std::pow(sum, scalar_type{1} / p);
  }

  template<typename Other>
    requires compatible_vectors_c_weak<raw_vector, Other>
  constexpr raw_vector cross(const Other& other) const
  { return lam::linalg::cross(*this, other); }

  constexpr raw_vector unit() const { return lam::linalg::unit(*this); }

  template<typename Other>
    requires compatible_vectors_c_weak<raw_vector, Other>
  constexpr bool operator==(const Other& other) const
  {
    if (m_dimension != other.size())
      return false;
    auto rhs = other.as_span();
    for (size_type i = 0; i < m_dimension; ++i)
    {
      if (m_data[i] != rhs[i])
        return false;
    }
    return true;
  }

  constexpr raw_vector operator-() const
  {
    raw_vector res(*this, m_alloc);
    for (size_type i = 0; i < m_dimension; ++i)
      res[i] = -res[i];
    return res;
  }

  template<typename Other>
    requires compatible_vectors_c_weak<raw_vector, Other>
  constexpr raw_vector& operator+=(const Other& other)
  {
    if (m_dimension != other.size())
      throw vector_exception::plus_equals_unequal_dim();
    auto rhs = other.as_span();
    for (size_type i = 0; i < m_dimension; ++i)
      m_data[i] += rhs[i];
    return *this;
  }

  constexpr pointer begin() noexcept { return m_data; }
  constexpr pointer end() noexcept { return m_data + m_dimension; }
  constexpr const_pointer begin() const noexcept { return m_data; }
  constexpr const_pointer end() const noexcept { return m_data + m_dimension; }
  constexpr const_pointer cbegin() const noexcept { return m_data; }
  constexpr const_pointer cend() const noexcept { return m_data + m_dimension; }

  constexpr std::span<scalar_type> as_span() noexcept { return {m_data, m_dimension}; }
  constexpr std::span<const scalar_type> as_span() const noexcept { return {m_data, m_dimension}; }

  template<typename Other>
    requires compatible_vectors_c_weak<raw_vector, Other>
  constexpr raw_vector& operator-=(const Other& other)
  {
    if (m_dimension != other.size())
      throw vector_exception::plus_equals_unequal_dim();
    auto rhs = other.as_span();
    for (size_type i = 0; i < m_dimension; ++i)
      m_data[i] -= rhs[i];
    return *this;
  }

  constexpr raw_vector& operator*=(const scalar_type& s)
  {
    for (size_type i = 0; i < m_dimension; ++i)
      m_data[i] *= s;
    return *this;
  }

  constexpr raw_vector& operator/=(const scalar_type& s)
  {
    if (s == scalar_type{0})
      throw vector_exception::div_by_zero();
    scalar_type inv_s = scalar_type{1} / s;
    for (size_type i = 0; i < m_dimension; ++i)
      m_data[i] *= inv_s;
    return *this;
  }

  constexpr raw_vector& normalize()
  {
    scalar_type n = norm();
    if (n == scalar_type{0})
      throw vector_exception::div_by_zero();
    scalar_type inv_n = scalar_type{1} / n;
    operator*=(inv_n); // Reuse operator
    return *this;
  }

  template<typename Other>
    requires compatible_vectors_c_weak<raw_vector, Other>
  constexpr scalar_type angle(const Other& other) const
  { return lam::linalg::angle(*this, other); }
  template<typename Other>
    requires compatible_vectors_c_weak<raw_vector, Other>
  constexpr raw_vector project(const Other& onto) const
  { return lam::linalg::project(*this, onto); }
  template<typename Other>
    requires compatible_vectors_c_weak<raw_vector, Other>
  constexpr raw_vector reject(const Other& from) const
  { return lam::linalg::reject(*this, from); }
  template<typename Other>
    requires compatible_vectors_c_weak<raw_vector, Other>
  constexpr raw_vector reflect(const Other& normal) const
  { return lam::linalg::reflect(*this, normal); }
  template<typename Other>
    requires compatible_vectors_c_weak<raw_vector, Other>
  constexpr raw_vector lerp(const Other& other, scalar_type t) const
  { return lam::linalg::lerp(*this, other, t); }
  template<typename Other>
    requires compatible_vectors_c_weak<raw_vector, Other>
  constexpr scalar_type distance(const Other& other) const
  { return lam::linalg::distance(*this, other); }
  template<typename Other>
    requires compatible_vectors_c_weak<raw_vector, Other>
  constexpr bool is_parallel(const Other& other, scalar_type tolerance = default_angle_tolerance<scalar_type>) const
  { return lam::linalg::is_parallel(*this, other, tolerance); }
  template<typename Other>
    requires compatible_vectors_c_weak<raw_vector, Other>
  constexpr bool is_orthogonal(const Other& other, scalar_type tolerance = default_angle_tolerance<scalar_type>) const
  { return lam::linalg::is_orthogonal(*this, other, tolerance); }
  template<typename B, typename C>
    requires compatible_vectors_c_weak<raw_vector, B> && compatible_vectors_c_weak<raw_vector, C>
  constexpr scalar_type triple_product(const B& b, const C& c) const
  { return lam::linalg::triple_product(*this, b, c); }
};

// Thin wrappers over the generic algorithms, fixed to `raw_vector`.

// Operands may carry different allocators; the result takes the LEFT one's.
export template<typename T, typename A1, typename A2>
constexpr raw_vector<T, A1> operator+(const raw_vector<T, A1>& a, const raw_vector<T, A2>& b)
{
  raw_vector<T, A1> res(a, a.get_allocator());
  res += b;
  return res;
}

// Operands may carry different allocators; the result takes the LEFT one's.
export template<typename T, typename A1, typename A2>
constexpr raw_vector<T, A1> operator-(const raw_vector<T, A1>& a, const raw_vector<T, A2>& b)
{
  raw_vector<T, A1> res(a, a.get_allocator());
  res -= b;
  return res;
}

export template<typename T, typename Alloc>
constexpr raw_vector<T, Alloc> operator*(const T& s, const raw_vector<T, Alloc>& v)
{
  raw_vector<T, Alloc> res(v, v.get_allocator());
  res *= s;
  return res;
}

export template<typename T, typename Alloc>
constexpr raw_vector<T, Alloc> operator*(const raw_vector<T, Alloc>& v, const T& s)
{ return s * v; }

export template<typename T, typename Alloc>
constexpr raw_vector<T, Alloc> operator/(const raw_vector<T, Alloc>& v, const T& s)
{
  raw_vector<T, Alloc> res(v, v.get_allocator());
  res /= s;
  return res;
}

// Spares the caller from naming an allocator type: a memory_resource* converts to the allocator.
export namespace pmr
{
template<typename T>
using raw_vector = lam::linalg::raw_vector<T, std::pmr::polymorphic_allocator<T>>;
} // namespace pmr

} // namespace lam::linalg

export template<typename T, typename Alloc>
struct std::formatter<lam::linalg::raw_vector<T, Alloc>>
{
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  auto format(const lam::linalg::raw_vector<T, Alloc>& v, std::format_context& ctx) const
  {
    auto out = ctx.out();
    std::format_to(out, "(");
    for (std::size_t i = 0; i < v.size(); ++i)
    {
      if (i > 0)
        std::format_to(out, ", ");
      std::format_to(out, "{}", v[i]);
    }
    return std::format_to(out, ")");
  }
};
