export module lam.linearalgebra:vectorspace.r_vector;

import std;
import lam.concepts;
import :vectorspace.concepts;
import :vectorspace.exceptions;
import :vectorspace.algorithms;

namespace lam::linalg
{

export 
template<typename T, typename Alloc = std::allocator<T>>
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
  struct Deleter
  {
    Alloc alloc;
    size_type n;
    constexpr Deleter(const Alloc& a = Alloc(), size_type s = 0) : alloc(a), n(s) {}
    constexpr void operator()(pointer p) 
    {
      if (p)
      {
        // Range-based destruction
        // std::destroy_n is standard range equivalent for raw pointers
        // But we need to use allocator traits destroy.
        // Range-based allocator destruction is not standard yet.
        // We will stick to ranges::for_each over a transformed view?
        // Or simpler: ranges::for_each on span.
        // Note: destroy expects valid pointer.
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

  constexpr explicit r_vector(size_type dim, const allocator_type& alloc = allocator_type())
    : m_alloc{alloc}, m_dimension{dim}, arrow{nullptr, Deleter{alloc, dim}}
  {
    if (dim == 0)
      throw vector_exception::non_pos();
    auto* ptr = std::allocator_traits<Alloc>::allocate(m_alloc, dim);
    arrow.reset(ptr);
    
    // Range-based construction:
    // We want to construct `dim` elements.
    // std::ranges::uninitialized_default_construct_n?? 
    // BUT we must use Allocator.
    // So iterating over the allocated range is best.
    std::ranges::for_each(std::span{ptr, dim}, [&](auto& e){
         std::allocator_traits<Alloc>::construct(m_alloc, &e);
    });
  }

  constexpr r_vector(std::initializer_list<scalar_type> il, const allocator_type& alloc = allocator_type())
    : m_alloc{alloc}, m_dimension{il.size()}, arrow{nullptr, Deleter{alloc, il.size()}}
  {
    if (il.size() == 0)
      throw vector_exception::non_pos();
    auto* ptr = std::allocator_traits<Alloc>::allocate(m_alloc, m_dimension);
    arrow.reset(ptr);
    
    // Zip init list and allocated memory
    // Note: C++23 allows zip view over init_list usually, or we just use iterators.
    // ranges::for_each(zip(il, span))
    for (auto&& [val, dest] : std::views::zip(il, std::span{ptr, m_dimension})) {
         std::allocator_traits<Alloc>::construct(m_alloc, &dest, val);
    }
  }

  constexpr r_vector(const std::vector<scalar_type>& v, const allocator_type& alloc = allocator_type())
    : m_alloc{alloc}, m_dimension{v.size()}, arrow{nullptr, Deleter{alloc, v.size()}}
  {
    if (v.size() == 0)
      throw vector_exception::non_pos();
    auto* ptr = std::allocator_traits<Alloc>::allocate(m_alloc, m_dimension);
    arrow.reset(ptr);
    
    for (auto&& [val, dest] : std::views::zip(v, std::span{ptr, m_dimension})) {
         std::allocator_traits<Alloc>::construct(m_alloc, &dest, val);
    }
  }

  r_vector(r_vector&&) noexcept = default;
  r_vector& operator=(r_vector&&) noexcept = default;

  constexpr r_vector(const r_vector& other) 
    : m_alloc{other.m_alloc},
      m_dimension{other.m_dimension}, arrow{nullptr, Deleter{m_alloc, other.m_dimension}}
  {
    if (m_dimension > 0)
    {
      auto* ptr = std::allocator_traits<Alloc>::allocate(m_alloc, m_dimension);
      arrow.reset(ptr);
      
      for (auto&& [val, dest] : std::views::zip(other.as_span(), std::span{ptr, m_dimension})) {
         std::allocator_traits<Alloc>::construct(m_alloc, &dest, val);
      }
    }
  }

  constexpr r_vector(const r_vector& other, const allocator_type& alloc)
    : m_alloc{alloc}, m_dimension{other.m_dimension}, arrow{nullptr, Deleter{alloc, other.m_dimension}}
  {
    if (m_dimension > 0)
    {
      auto* ptr = std::allocator_traits<Alloc>::allocate(m_alloc, m_dimension);
      arrow.reset(ptr);
      
      for (auto&& [val, dest] : std::views::zip(other.as_span(), std::span{ptr, m_dimension})) {
         std::allocator_traits<Alloc>::construct(m_alloc, &dest, val);
      }
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

  // Generic Delegates
  constexpr scalar_type dot(const r_vector& other) const { return lam::linalg::dot(*this, other); }
  constexpr scalar_type norm2() const { return lam::linalg::norm2(*this); }
  constexpr scalar_type norm() const { return lam::linalg::norm(*this); }

  constexpr scalar_type pnorm(scalar_type p) const
  {
    // Ranges pnorm!
    // sum = fold_left( transform(abs(v), pow(p)), 0, + )
    // Not standard, but doable.
    // C++23 fold_left requires <algorithm>
    auto abs_pow = std::views::transform(as_span(), [p](auto val){ return std::pow(std::abs(val), p); });
    scalar_type sum = std::ranges::fold_left(abs_pow, scalar_type{0}, std::plus<scalar_type>{});
    return std::pow(sum, scalar_type{1} / p);
  }

  constexpr r_vector cross(const r_vector& other) const { return lam::linalg::cross(*this, other); }
  constexpr r_vector unit() const { return lam::linalg::unit(*this); }

  constexpr bool operator==(const r_vector& other) const
  {
    if (m_dimension != other.m_dimension) return false;
    return std::ranges::equal(as_span(), other.as_span());
  }

  constexpr r_vector operator-() const
  {
    r_vector res(*this, m_alloc);
    // Range transform in place?
    std::ranges::transform(res.as_span(), res.begin(), [](auto val){ return -val; });
    return res;
  }

  constexpr r_vector& operator+=(const r_vector& other)
  {
    if (m_dimension != other.m_dimension)
      throw vector_exception::plus_equals_unequal_dim();
    
    // Zip + for_each to add
    for (auto&& [a, b] : std::views::zip(as_span(), other.as_span())) {
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

  constexpr r_vector& operator-=(const r_vector& other)
  {
    if (m_dimension != other.m_dimension)
      throw vector_exception::plus_equals_unequal_dim();
    
    for (auto&& [a, b] : std::views::zip(as_span(), other.as_span())) {
        a -= b;
    }
    return *this;
  }
  
  constexpr r_vector& operator*=(const scalar_type& s)
  {
    // transform in place
    // Or for_each
    std::ranges::for_each(as_span(), [s](auto& val){ val *= s; });
    return *this;
  }

  constexpr r_vector& operator/=(const scalar_type& s)
  {
    if (s == scalar_type{0})
      throw vector_exception::div_by_zero();
    scalar_type inv_s = scalar_type{1} / s;
    std::ranges::for_each(as_span(), [inv_s](auto& val){ val *= inv_s; });
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
  constexpr scalar_type angle(const r_vector& other) const { return lam::linalg::angle(*this, other); }
  constexpr r_vector project(const r_vector& onto) const { return lam::linalg::project(*this, onto); }
  constexpr r_vector reject(const r_vector& from) const { return lam::linalg::reject(*this, from); }
  constexpr r_vector reflect(const r_vector& normal) const { return lam::linalg::reflect(*this, normal); }
  constexpr r_vector lerp(const r_vector& other, scalar_type t) const { return lam::linalg::lerp(*this, other, t); }
  constexpr scalar_type distance(const r_vector& other) const { return lam::linalg::distance(*this, other); }
  constexpr bool is_parallel(const r_vector& other, scalar_type tolerance = scalar_type{1e-10}) const { return lam::linalg::is_parallel(*this, other, tolerance); }
  constexpr bool is_orthogonal(const r_vector& other, scalar_type tolerance = scalar_type{1e-10}) const { return lam::linalg::is_orthogonal(*this, other, tolerance); }
  constexpr scalar_type triple_product(const r_vector& b, const r_vector& c) const { return lam::linalg::triple_product(*this, b, c); }
};

// Operators
export template<typename T, typename Alloc>
constexpr r_vector<T, Alloc> operator+(const r_vector<T, Alloc>& a, const r_vector<T, Alloc>& b)
{
  r_vector<T, Alloc> res(a, a.get_allocator());
  res += b;
  return res;
}

export template<typename T, typename Alloc>
constexpr r_vector<T, Alloc> operator-(const r_vector<T, Alloc>& a, const r_vector<T, Alloc>& b)
{
  r_vector<T, Alloc> res(a, a.get_allocator());
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
{
  return s * v;
}

export template<typename T, typename Alloc>
constexpr r_vector<T, Alloc> operator/(const r_vector<T, Alloc>& v, const T& s)
{
  r_vector<T, Alloc> res(v, v.get_allocator());
  res /= s;
  return res;
}

} // namespace lam::linalg
