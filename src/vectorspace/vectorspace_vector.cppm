export module lam.linearalgebra:vectorspace.vector;

import std;
import lam.concepts;
import :vectorspace.concepts;
import :vectorspace.exceptions;
import :vectorspace.algorithms;

namespace lam::linalg
{

export template<typename T, typename Alloc = std::allocator<T>>
  requires lam::concepts::experimental::ring_element_c_weak<T>
class vector
{
public:
  using scalar_type = T;
  using allocator_type = Alloc;
  using size_type = std::size_t;
  using pointer = scalar_type*;
  using const_pointer = const scalar_type*;

private:
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
          std::allocator_traits<Alloc>::destroy(alloc, p + i);
        std::allocator_traits<Alloc>::deallocate(alloc, p, n);
      }
    }
  };

  [[no_unique_address]] Alloc m_alloc;
  size_type m_dimension;
  std::unique_ptr<scalar_type[], Deleter> arrow;

public:
  constexpr vector() : m_alloc{}, m_dimension{0}, arrow{nullptr, Deleter{}} {}

  constexpr explicit vector(size_type dim, const allocator_type& alloc = allocator_type())
    : m_alloc{alloc}, m_dimension{dim}, arrow{nullptr, Deleter{alloc, dim}}
  {
    if (dim == 0)
      throw vector_exception::non_pos();
    auto* ptr = std::allocator_traits<Alloc>::allocate(m_alloc, dim);
    arrow.reset(ptr);
    for (size_type i = 0; i < dim; ++i)
    {
      std::allocator_traits<Alloc>::construct(m_alloc, ptr + i);
    }
  }

  constexpr vector(std::initializer_list<scalar_type> il, const allocator_type& alloc = allocator_type())
    : m_alloc{alloc}, m_dimension{il.size()}, arrow{nullptr, Deleter{alloc, il.size()}}
  {
    if (il.size() == 0)
      throw vector_exception::non_pos();
    auto* ptr = std::allocator_traits<Alloc>::allocate(m_alloc, m_dimension);
    arrow.reset(ptr);
    size_type i = 0;
    for (const auto& val : il)
    {
      std::allocator_traits<Alloc>::construct(m_alloc, ptr + i, val);
      ++i;
    }
  }

  constexpr vector(const std::vector<scalar_type>& v, const allocator_type& alloc = allocator_type())
    : m_alloc{alloc}, m_dimension{v.size()}, arrow{nullptr, Deleter{alloc, v.size()}}
  {
    if (v.size() == 0)
      throw vector_exception::non_pos();
    auto* ptr = std::allocator_traits<Alloc>::allocate(m_alloc, m_dimension);
    arrow.reset(ptr);
    for (size_type i = 0; i < m_dimension; ++i)
    {
      std::allocator_traits<Alloc>::construct(m_alloc, ptr + i, v[i]);
    }
  }

  vector(vector&&) noexcept = default;
  vector& operator=(vector&&) noexcept = default;

  constexpr vector(const vector& other)
    : m_alloc{other.m_alloc}, m_dimension{other.m_dimension}, arrow{nullptr, Deleter{m_alloc, other.m_dimension}}
  {
    if (m_dimension > 0)
    {
      auto* ptr = std::allocator_traits<Alloc>::allocate(m_alloc, m_dimension);
      arrow.reset(ptr);
      for (size_type i = 0; i < m_dimension; ++i)
      {
        std::allocator_traits<Alloc>::construct(m_alloc, ptr + i, other[i]);
      }
    }
  }

  constexpr vector(const vector& other, const allocator_type& alloc)
    : m_alloc{alloc}, m_dimension{other.m_dimension}, arrow{nullptr, Deleter{alloc, other.m_dimension}}
  {
    if (m_dimension > 0)
    {
      auto* ptr = std::allocator_traits<Alloc>::allocate(m_alloc, m_dimension);
      arrow.reset(ptr);
      for (size_type i = 0; i < m_dimension; ++i)
      {
        std::allocator_traits<Alloc>::construct(m_alloc, ptr + i, other[i]);
      }
    }
  }

  constexpr vector& operator=(const vector& other)
  {
    if (this != &other)
    {
      vector tmp(other, m_alloc);
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

  constexpr scalar_type dot(const vector& other) const { return lam::linalg::dot(*this, other); }
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

  constexpr vector cross(const vector& other) const { return lam::linalg::cross(*this, other); }

  constexpr vector unit() const { return lam::linalg::unit(*this); }

  constexpr bool operator==(const vector& other) const
  {
    if (m_dimension != other.m_dimension)
      return false;
    for (size_type i = 0; i < m_dimension; ++i)
    {
      if ((*this)[i] != other[i])
        return false;
    }
    return true;
  }

  constexpr vector operator-() const
  {
    vector res(*this, m_alloc);
    for (size_type i = 0; i < m_dimension; ++i)
      res[i] = -res[i];
    return res;
  }

  constexpr vector& operator+=(const vector& other)
  {
    if (m_dimension != other.m_dimension)
      throw vector_exception::plus_equals_unequal_dim();
    for (size_type i = 0; i < m_dimension; ++i)
      arrow[i] += other.arrow[i];
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

  constexpr vector& operator-=(const vector& other)
  {
    if (m_dimension != other.m_dimension)
      throw vector_exception::plus_equals_unequal_dim();
    for (size_type i = 0; i < m_dimension; ++i)
      arrow[i] -= other.arrow[i];
    return *this;
  }

  constexpr vector& operator*=(const scalar_type& s)
  {
    for (size_type i = 0; i < m_dimension; ++i)
      arrow[i] *= s;
    return *this;
  }

  constexpr vector& operator/=(const scalar_type& s)
  {
    if (s == scalar_type{0})
      throw vector_exception::div_by_zero();
    scalar_type inv_s = scalar_type{1} / s;
    for (size_type i = 0; i < m_dimension; ++i)
      arrow[i] *= inv_s;
    return *this;
  }

  constexpr vector& normalize()
  {
    // Could define a generic normalize(v) that modifies in place?
    // Current generic algorithms return new objects (functional).
    // Let's stick to implementation here or define 'normalize_inplace' generic.
    // For now, keep local logic or use: *this = unit(); (expensive copy)
    // Optimization: logic is simple enough to keep or we add generic 'normalize_inplace'.
    // Sticking to local optimization to avoid copy assignment if possible.
    scalar_type n = norm();
    if (n == scalar_type{0})
      throw vector_exception::div_by_zero();
    scalar_type inv_n = scalar_type{1} / n;
    operator*=(inv_n); // Reuse operator
    return *this;
  }

  // Delegate to generic algorithms
  constexpr scalar_type angle(const vector& other) const { return lam::linalg::angle(*this, other); }
  constexpr vector project(const vector& onto) const { return lam::linalg::project(*this, onto); }
  constexpr vector reject(const vector& from) const { return lam::linalg::reject(*this, from); }
  constexpr vector reflect(const vector& normal) const { return lam::linalg::reflect(*this, normal); }
  constexpr vector lerp(const vector& other, scalar_type t) const { return lam::linalg::lerp(*this, other, t); }
  constexpr scalar_type distance(const vector& other) const { return lam::linalg::distance(*this, other); }
  constexpr bool is_parallel(const vector& other, scalar_type tolerance = scalar_type{1e-10}) const
  {
    return lam::linalg::is_parallel(*this, other, tolerance);
  }
  constexpr bool is_orthogonal(const vector& other, scalar_type tolerance = scalar_type{1e-10}) const
  {
    return lam::linalg::is_orthogonal(*this, other, tolerance);
  }
  constexpr scalar_type triple_product(const vector& b, const vector& c) const
  {
    return lam::linalg::triple_product(*this, b, c);
  }
};

// Operators remain as thin wrappers or can also be genericized?
// For THIS file (vector class partition), they operate on `vector`.

export template<typename T, typename Alloc>
constexpr vector<T, Alloc> operator+(const vector<T, Alloc>& a, const vector<T, Alloc>& b)
{
  vector<T, Alloc> res(a, a.get_allocator());
  res += b;
  return res;
}

export template<typename T, typename Alloc>
constexpr vector<T, Alloc> operator-(const vector<T, Alloc>& a, const vector<T, Alloc>& b)
{
  vector<T, Alloc> res(a, a.get_allocator());
  res -= b;
  return res;
}

export template<typename T, typename Alloc>
constexpr vector<T, Alloc> operator*(const T& s, const vector<T, Alloc>& v)
{
  vector<T, Alloc> res(v, v.get_allocator());
  res *= s;
  return res;
}

export template<typename T, typename Alloc>
constexpr vector<T, Alloc> operator*(const vector<T, Alloc>& v, const T& s)
{
  return s * v;
}

export template<typename T, typename Alloc>
constexpr vector<T, Alloc> operator/(const vector<T, Alloc>& v, const T& s)
{
  vector<T, Alloc> res(v, v.get_allocator());
  res /= s;
  return res;
}

// Re-export specific generic algorithms if needed, OR relies on users importing :algorithms?
// The primary `vectorspace` module will export algorithms.
// We keep free function wrappers here for ADL if desired?
// No, removing specific wrappers (dot, cross, etc) from here because they are now in algorithms partition
// AND exported by the primary module.
// BUT `vector` users expect `dot(v1, v2)` to work via ADL?
// The algorithms in `lam::linalg` namespace will be found.

} // namespace lam::linalg

export template<typename T, typename Alloc>
struct std::formatter<lam::linalg::vector<T, Alloc>>
{
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  auto format(const lam::linalg::vector<T, Alloc>& v, std::format_context& ctx) const
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
