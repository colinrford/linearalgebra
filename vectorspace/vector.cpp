
#include "vector.hpp"
#include <algorithm>
#include <cstddef>
#include <chrono>
#include <compare>
#include <execution>
#include <iterator>
#include <ranges>
#include <utility>

namespace linalg
{

/*auto make_indexing_set = [](const std::size_t n) -> std::list<std::size_t>
{
	std::list<std::size_t> ell(n);
	std::iota(ell.begin(), ell.end(), 0);
	return ell;
}; */

vector::vector(const std::size_t dim)
      : dimension{dim}
{
  if (dim != 0) [[likely]]
  {
    arrow = std::make_unique<double[]>(dimension);
  } else  throw vector_exception::non_pos();
}

// dangerous, now a private constructor
vector::vector(const std::size_t dim, std::unique_ptr<double[]> elem)
      : dimension{dim}
{
  if (elem == nullptr) [[unlikely]]
    throw vector_exception::null_construction();
  if (dim == 0) [[unlikely]]
    throw vector_exception::non_pos();
  arrow = std::move(elem);
}

vector::vector(const std::vector<double> elem)
      : dimension{elem.size()}
{
  if (dimension > 0) [[likely]]
  {
    arrow = std::make_unique<double[]>(dimension);
    for (std::size_t i = 0; i < dimension; i++)
      arrow[i] = elem[i];
  } else  throw vector_exception::non_pos();
}

vector::vector(const std::initializer_list<double> il)
      : dimension{il.size()}
{
  if (il.size() == 0) [[unlikely]]
    throw vector_exception::non_pos();
  arrow = std::make_unique<double[]>(il.size());
  auto arrow_ptr = arrow.get();
  for (auto entry : il)
  {
    *arrow_ptr = entry;
    arrow_ptr++;
  }
}

vector::vector(vector&& v) noexcept
      : dimension{v.dimension},
        arrow{std::move_if_noexcept(v.arrow)}
{
  v.dimension = 0;
}

vector& vector::operator=(vector&& v) noexcept
{
  this->arrow = std::move(v.arrow);
  this->dimension = v.dimension;
  v.dimension = 0;
  return *this;
}

auto vector::operator[](const std::size_t index)
                      -> decltype(arrow[index])
{
  if (index < this->dimension)
    return arrow[index];
  else
    throw vector_exception::out_of_bounds();
}

auto vector::operator[](const std::size_t index) const
                      -> decltype(arrow[index])
{
  if (index < this->dimension)
    return arrow[index];
  else
    throw vector_exception::out_of_bounds();
}

vector::~vector() = default;

/*
// not sure why i have these besides having triple comparison for the hekovit
auto vector::operator[](const int index)
                      -> decltype(arrow[index])
{
  if (0 <= index < this->dimension)
    return arrow[index];
  else
    throw vector_exception::out_of_bounds();
}
// not sure why i have these besides having triple comparison for the hekovit
auto vector::operator[](const int index) const
                      -> decltype(arrow[index])
{
  if (0 <= index < this->dimension)
    return arrow[index];
  else
    throw vector_exception::out_of_bounds();
} */

/*
 *  iterator implementation
 */
struct vector::iterator
{
public:

  using value_type = double;
  using difference_type = std::ptrdiff_t;
  using pointer = double*;
  using reference = double&;
  using iterator_category = std::random_access_iterator_tag;

public:

  constexpr iterator() noexcept = default;

  constexpr iterator(pointer arrow_ptr) noexcept { iter_ptr = arrow_ptr; }

  constexpr pointer operator->() noexcept { return iter_ptr; }

  constexpr reference operator*() noexcept { return *iter_ptr; }

  constexpr reference operator[](difference_type index)
  { return *(iter_ptr + index); }

  constexpr iterator& operator++()
  {
    ++iter_ptr;
    return *this;
  }

  constexpr iterator operator++(int)
  {
    auto iter = *this;
    ++(*this);
    return iter;
  }

  constexpr iterator& operator--()
  {
    --iter_ptr;
    return *this;
  }

  constexpr iterator operator--(int)
  {
    auto iter = *this;
    --(*this);
    return iter;
  }

  constexpr iterator& operator+=(int off)
  {
    iter_ptr += off;
    return *this;
  }

  constexpr iterator operator+(difference_type& off)
  {
    auto iter = *this;
    return iter += off;// return iterator(iter + off);// based off gcc stl impl
  }

  constexpr iterator& operator-=(int off) { return *this += -off; }

  constexpr iterator operator-(difference_type& off)
  {
    auto iter = *this;
    return iter -= off;
  }

  constexpr bool operator==(const iterator& other) const noexcept
  { return iter_ptr == other.iter_ptr; }

  constexpr bool operator!=(const iterator& other) const noexcept
  { return iter_ptr != other.iter_ptr; } // should it be !( == )?

  constexpr bool operator<(iterator other) const noexcept
  { return iter_ptr < other.iter_ptr; }

  constexpr bool operator>(iterator other) const noexcept
  { return iter_ptr > other.iter_ptr; }

  constexpr bool operator<=(iterator other) const noexcept
  { return iter_ptr <= other.iter_ptr; }

  constexpr bool operator>=(iterator other) const noexcept
  { return iter_ptr >= other.iter_ptr; }

  constexpr difference_type operator-(iterator other) const noexcept
  { return iter_ptr - other.iter_ptr; }

  constexpr friend auto operator <=>(iterator, iterator) = default;

  /*\constexpr friend difference_type operator-(iterator first, iterator second)
  { return *first - *second; }*/ // are

  constexpr friend iterator operator+(iterator it, int off)
  { return it += off; } // these

  constexpr friend iterator operator-(iterator it, int off)
  { return it -= off; } // okay? sure

  //friend iterator operator+(int off, iterator);

private:
  pointer iter_ptr;
};

vector::iterator vector::begin() { return vector::iterator(arrow.get()); }

vector::iterator vector::end()
{ return vector::iterator(arrow.get() + dimension); }

vector::iterator begin(vector& v) { return v.begin(); }

vector::iterator end(vector& v) { return v.end(); }
/*
 *    END iterator  *
 */

/*
*    BEGIN const_iterator    *
*/
struct vector::const_iterator
{
public:

 using value_type = const double;
 using difference_type = std::ptrdiff_t;
 using pointer = const double*;
 using reference = const double&;
 using iterator_category = std::random_access_iterator_tag;

public:

 constexpr const_iterator() noexcept = default;

 constexpr const_iterator(pointer arrow_ptr) noexcept
 { iter_ptr = arrow_ptr; }

 constexpr pointer operator->() const noexcept { return iter_ptr; }

 constexpr reference operator*() const noexcept { return *iter_ptr; }

 constexpr reference operator[](difference_type index) const
 { return *(iter_ptr + index); }

 constexpr const_iterator& operator++()
 {
   ++iter_ptr;
   return *this;
 }

 constexpr const_iterator operator++(int)
 {
   auto iter = *this;
   ++(*this);
   return iter;
 }

 constexpr const_iterator& operator--()
 {
   --iter_ptr;
   return *this;
 }

 constexpr const_iterator operator--(int)
 {
   auto iter = *this;
   --(*this);
   return iter;
 }

 constexpr const_iterator& operator+=(const int off)
 {
   iter_ptr += off;
   return *this;
 }

 constexpr const_iterator operator+(difference_type& off)
 {
   auto iter = *this;
   return iter += off;// return iterator(iter + off);// based off gcc stl impl
 }

 constexpr const_iterator& operator-=(int off) { return *this += -off; }

 constexpr const_iterator operator-(difference_type& off)
 {
   auto iter = *this;
   return iter -= off;
 }

 constexpr bool operator==(const const_iterator& other) const noexcept
 { return iter_ptr == other.iter_ptr; }

 constexpr bool operator!=(const const_iterator& other) const noexcept
 { return iter_ptr != other.iter_ptr; }

 constexpr bool operator<(const_iterator other) const noexcept
 { return iter_ptr < other.iter_ptr; }

 constexpr bool operator>(const_iterator other) const noexcept
 { return iter_ptr > other.iter_ptr; }

 constexpr bool operator<=(const_iterator other) const noexcept
 { return iter_ptr <= other.iter_ptr; }

 constexpr bool operator>=(const_iterator other) const noexcept
 { return iter_ptr >= other.iter_ptr; }

 /*constexpr difference_type operator-(const_iterator other) const noexcept
 { return iter_ptr - other.iter_ptr; } */

 constexpr friend auto operator <=>(const_iterator, const_iterator) = default;

 constexpr friend double operator-(const_iterator first,
                                   const_iterator second)
 { return *first - *second; }

 constexpr friend const_iterator operator+(const_iterator it, int off)
 { return it += off; }

 constexpr friend const_iterator operator-(const_iterator it, int off)
 { return it -= off; }

private:
 pointer iter_ptr;
};

vector::const_iterator vector::begin() const
{ return vector::const_iterator(arrow.get()); }

vector::const_iterator vector::end() const
{ return vector::const_iterator(arrow.get() + dimension); }

const vector::const_iterator vector::cbegin() const
{ return vector::const_iterator(arrow.get()); }

const vector::const_iterator vector::cend() const
{ return vector::const_iterator(arrow.get() + dimension); }

/*
*  END const_iterator implementation
*/

void vector::give_label(std::string s)
{ label = std::make_optional(s); }

void vector::give_label(std::string_view s)
{ label = std::make_optional(s); }

std::string_view vector::get_label() const
{
  if (label.has_value())
    return *label;
  else
    return "";
}

std::string_view vector::name() const
{ return this->get_label(); }

void vector::compute_info()
{
  if (!info.has_value())
  {
    const auto scoped_norm_sq = this->dot(*this);
    const auto scoped_norm = sqrt(scoped_norm_sq);
    const auto scoped_info = vector_info(scoped_norm, scoped_norm_sq);
    auto scoped_pointed_info = std::make_unique<vector_info>(scoped_info);
    info = std::make_optional(std::move(scoped_pointed_info));
  }
}

// should these be further restricted or removed (at least as public members?)
double* vector::get() const
{ return arrow.get(); }

double* vector::get_arrow() const
{ return arrow.get(); }

double vector::dot(const vector& v2) const
{
  const auto d1 = this->dimension;
  const auto d2 = v2.dimension;
  std::size_t min = 0;
  if (d1 != d2)
    min = std::min(d1, d2);
  else
    min = d1;

  double dot = 0.0;
  for (std::size_t i = 0; i < min; i++)
    dot += this->arrow[i] * v2.arrow[i];

  return (!compare(dot, 0.)) ? dot : 0.;
}

double vector::better_dot(const vector& v2) const
{
  double dot = std::inner_product(this->begin(), this->end(),
                                  v2.begin(), 0.);
  // this is actually not necessarily faster (at least, without
  // an optiminization flag) general in my limited test runs
  // for large n: via gcc it's faster; via clang it's slower
  return (!compare(dot, 0.)) ? dot : 0.;
}

double vector::even_better_dot(const vector& v2) const
{
  double dot = std::transform_reduce(this->begin(), this->end(),
                                     v2.begin(), 0.);
  // read previous function
  return (!compare(dot, 0.)) ? dot : 0.;
}

double dot(const vector& v1, const vector& v2)
{ return v1.dot(v2); }

double vector::mag()
{
  if (!info.has_value())
    this->compute_info();
  return (*info)->norm;
}

double vector::mag() const
{
  if (info.has_value())
    return (*info)->norm;
  else
    return sqrt(this->dot(*this));
}

double vector::mag2()
{
  if (!info.has_value())
    this->compute_info();
  return (*info)->norm_squared;
}

double vector::mag2() const
{
  if (info.has_value())
    return (*info)->norm_squared;
  else
    return this->dot(*this);
}

double vector::norm()
{
  if (!info.has_value())
    this->compute_info();
  return (*info)->norm;
}

double vector::norm() const
{
  if (info.has_value())
    return (*info)->norm;
  else
    return sqrt(this->dot(*this));
}

double vector::norm2()
{
  if (!info.has_value())
    this->compute_info();
  return (*info)->norm_squared;
}

double vector::norm2() const
{
  if (info.has_value())
    return (*info)->norm_squared;
  else
    return this->dot(*this);
}

vector vector::cross(const vector& v2) const
{
  const auto d1 = this->dimension;
  const auto d2 = v2.dimension;
  if ((d1 != d2) || (d1 != 3))
    throw vector_exception::cross_undef();

  auto elem = std::make_unique<double[]>(d1);
  elem[0] = this->arrow[1] * v2.arrow[2] - this->arrow[2] * v2.arrow[1];
  elem[1] = -(this->arrow[0] * v2.arrow[2] - this->arrow[2] * v2.arrow[0]);
  elem[2] = this->arrow[0] * v2.arrow[1] - this->arrow[1] * v2.arrow[0];

  for (std::size_t i = 0; i < d1; i++)
    if (compare(elem[i], 0.)) // expensive and possibly unnecessary?
      elem[i] = 0.;

  return vector(d1, std::move(elem));
}

vector cross(const vector& v1, const vector& v2)
{ return v1.cross(v2); }

vector vector::add(const vector& v2) const
{
  const auto d1 = this->dimension;
  const auto d2 = v2.dimension;
  if (d1 == d2)
  {
    auto elem = std::make_unique<double[]>(d1);
    for (std::size_t i = 0; i < d1; i++)
      elem[i] = this->arrow[i] + v2.arrow[i];
    return vector(d1, std::move(elem));
  } else
  {
    const auto [min, max] = std::minmax(d1, d2);
    auto elem = std::make_unique<double[]>(max);
    for (std::size_t i = 0; i < min; i++)
      elem[i] = this->arrow[i] + v2.arrow[i];
    if (d1 == min)
      for (std::size_t i = min; i < max; i++)
        elem[i] = v2.arrow[i];
    else
      for (std::size_t i = min; i < max; i++)
        elem[i] = this->arrow[i];
    return vector(max, std::move(elem));
  }
}

vector add(const vector& v1, const vector& v2)
{ return v1.add(v2); }

vector vector::subtract(const vector& v2) const
{
  const auto d1 = this->dimension;
  const auto d2 = v2.dimension;
  if (d1 == d2)
  {
    auto elem = std::make_unique<double[]>(d1);
    for (std::size_t i = 0; i < d1; i++)
      elem[i] = this->arrow[i] - v2.arrow[i];
    return vector(d1, std::move(elem));
  } else
  {
    const auto [min, max] = std::minmax(d1, d2);
    auto elem = std::make_unique<double[]>(max);
    for (std::size_t i = 0; i < min; i++)
      elem[i] = this->arrow[i] - v2.arrow[i];
    if (d1 == min)
      for (std::size_t i = min; i < max; i++)
        elem[i] = -v2.arrow[i];
    else
      for (std::size_t i = min; i < max; i++)
        elem[i] = this->arrow[i];
    return vector(max, std::move(elem));
  }
}

vector subtract(const vector& v1, const vector& v2)
{ return v1.subtract(v2); }

vector vector::scalar(const double s) const
{
  const std::size_t d = this->dimension;
  auto elem = std::make_unique<double[]>(d);
  for (std::size_t i = 0; i < d; i++)
    elem[i] = s * arrow[i];

  return vector(d, std::move(elem));
}

vector vector::unit() const
{
  double norm;
  if (info.has_value())
    norm = (*info)->norm;
  else
    norm = this->norm();

  if (compare(norm, 1.))
    return 1 * *this; // lazy, probably expensive copy

  const double one_over_norm = 1 / norm;
  const std::size_t d = this->dimension;
  auto elem = std::make_unique<double[]>(d);
  for (std::size_t i = 0; i < d; i++)
    elem[i] = this->arrow[i] * one_over_norm;

  return vector(d, std::move(elem));
}

constexpr bool vector::equals(const vector& v2) const
{
  const auto d1 = this->dimension;
  const auto d2 = v2.dimension;

  if (d1 != d2)
    return false;
  for (std::size_t i = 0; i < d1; i++)
    if (!compare(this->arrow[i], v2.arrow[i]))
      return false;

  return true;
}

vector& vector::operator+=(const vector& v2)
{
  const auto d1 = this->dimension;
  const auto d2 = v2.dimension;
  if (d1 == d2)
  {
    for (std::size_t i = 0; i < d1; i++)
    {
      std::cout << "v2.arrow[" << i << "] = " << v2.arrow[i] << std::endl;
      this->arrow[i] += v2.arrow[i];
    }
    return *this;
  } else  throw vector_exception::plus_equals_unequal_dim();
}

vector operator+(const vector& v1, const vector& v2) { return v1.add(v2); }

//vector operator+(const vector v1, const vector v2) { return v1.add(v2); }

vector operator-(const vector& v1, const vector& v2)
{ return v1.subtract(v2); }

vector operator*(const double d, const vector& v) { return v.scalar(d); }

vector operator*(const vector& v, const double d) { return v.scalar(d); }

vector operator/(const vector& v, const double d)
{
  if (compare(d, 0.))
    throw vector_exception::div_by_zero();
  return v.scalar(1 / d);
}

bool operator==(const vector& v1, const vector& v2)
{ return v1.equals(v2); }

/* bool compare(const double a, const double b)
{
  constexpr double epsilon = std::numeric_limits<double>::epsilon();

  if (abs(b - a) < epsilon)
    return true;
  else
    return false;
} */


double vector::pnorm(const double p) const
{
 //double little_lpnorm = 0.;
 const double one_over_p = 1 / p;
 /* for (const auto& component : *this)
 { little_lpnorm += std::pow(component, p); } */
 double little_lp_norm = std::accumulate(this->begin(), this->end(), 0.,
                                        [&] (auto powered, auto powerless) {
                                          return std::forward<double>(powered)
                                               + std::pow(std::abs(powerless),
                                                          p);
                                        });
 return std::pow(little_lp_norm, one_over_p);
}

double vector::lpnorm(const double p) const
{ return this->pnorm(p); }

vector make_zero_vector(const std::size_t dim)
{ return vector(dim); }

void vector::print() const
{
  std::cout << "(";
  auto i = this->begin();
  for (; i != this->end() - 1; ++i)
    std::cout << *i << ", ";
  std::cout << *i << ")" << std::endl;
}

void print(const vector& v)
{
  std::cout << "(";
  const auto all_but_last = std::ranges::subrange(v.begin(), v.end() - 1);
  for (const auto entry : all_but_last)
    std::cout << entry << ", ";
  std::cout << *(v.end() - 1) << ")" << std::endl;
}

}
