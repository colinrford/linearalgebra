
#include "vector.hpp"
#include <algorithm>
#include <cstddef>
#include <chrono>
#include <compare>
#include <iterator>
#include <utility>

namespace linalg
{

  auto make_indexing_set = [](const std::size_t n) -> std::list<std::size_t> {
  	std::list<std::size_t> ell(n);
  	std::iota(ell.begin(), ell.end(), 0);
  	return ell;
  }; //this may actually be too slow to use, sadly

vector::vector(const std::size_t dim)
{
  if (dim != 0)
  {
    dimension = dim;
    arrow = std::make_unique<double[]>(dimension);
    for (auto i = 0; i < dimension; i++)
      if (i != 0)
        arrow[i] = 0.0;
      else
        arrow[i] = 1.0;
  } else  throw vector_exception::non_pos();
}

vector::vector(const std::size_t dim, std::unique_ptr<double[]> elem)
{
  if (elem == nullptr)
    throw vector_exception::null_construction();
  if (dim == 0)
    throw vector_exception::non_pos();
  dimension = dim;
  arrow = std::move(elem);
}

vector::vector(const std::vector<double> elem)
{
  dimension = elem.size();
  if (dimension > 0)
  {
    arrow = std::make_unique<double[]>(dimension);
    for (auto i = 0; i < dimension; i++)
      arrow[i] = elem[i];
  } else  throw vector_exception::non_pos();
}

vector::vector(vector&& v) noexcept
      : arrow{std::move_if_noexcept(v.arrow)}, dimension{v.dimension}
{
  v.dimension = 0;
}

vector& vector::operator=(vector&& v)
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

constexpr std::size_t vector::get_dimension() const { return dimension; }

constexpr std::size_t vector::size() const { return dimension; }

constexpr std::size_t vector::length() const { return dimension; }

// Raw pointers are okay if they are non-owning
double* vector::get_arrow() const
{
  return arrow.get();
}

double vector::mag() const
{
  return sqrt(this->dot(*this));
}

double vector::mag2() const
{
  return this->dot(*this);
}

double vector::norm() const
{
  return sqrt(this->dot(*this));
}

double vector::norm2() const
{
  return this->dot(*this);
}

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
  for (auto i = 0; i < min; i++)
    dot += this->arrow[i] * v2.arrow[i];

  /*if (compare(dot, 0))
    dot = 0;*/

  return (!compare(dot, 0.)) ? dot : 0.;
}

double dot(vector& v1, const vector& v2)
{
  return v1.dot(v2);
}

vector vector::cross(const vector& v2)
{
  const auto d1 = this->dimension;
  const auto d2 = v2.dimension;
  if ((d1 != d2) || (d1 != 3))
    throw vector_exception::cross_undef();

  auto elem = std::make_unique<double[]>(d1);
  elem[0] = this->arrow[1] * v2.arrow[2] - this->arrow[2] * v2.arrow[1];
  elem[1] = -(this->arrow[0] * v2.arrow[2] - this->arrow[2] * v2.arrow[0]);
  elem[2] = this->arrow[0] * v2.arrow[1] - this->arrow[1] * v2.arrow[0];

  for (auto i = 0; i < d1; i++)
    if (compare(elem[i], 0.)) // expensive and possibly unnecessary?
      elem[i] = 0.;

  return vector(d1, std::move(elem));
}

vector cross(vector& v1, const vector& v2)
{
  return v1.cross(v2);
}

vector vector::add(const vector& v2)
{
  const auto d1 = this->dimension;
  const auto d2 = v2.dimension;
  if (d1 == d2)
  {
    auto elem = std::make_unique<double[]>(d1);
    for (auto i = 0; i < d1; i++)
      elem[i] = this->arrow[i] + v2.arrow[i];
    return vector(d1, std::move(elem));
  } else
  {
    const auto [min, max] = std::minmax(d1, d2);
    auto elem = std::make_unique<double[]>(max);
    for (auto i = 0; i < min; i++)
    {
      elem[i] = this->arrow[i] + v2.arrow[i];
    }
    if (d1 == min)
      for (auto i = min; i < max; i++)
        elem[i] = v2.arrow[i];
    else
      for (auto i = min; i < max; i++)
        elem[i] = this->arrow[i];
    return vector(max, std::move(elem));
  }
}

vector add(vector& v1, const vector& v2)
{
  return v1.add(v2);
}

vector vector::subtract(const vector& v2)
{
  const auto d1 = this->dimension;
  const auto d2 = v2.dimension;
  if (d1 == d2)
  {
    auto elem = std::make_unique<double[]>(d1);
    for (auto i = 0; i < d1; i++)
      elem[i] = this->arrow[i] - v2.arrow[i];
    return vector(d1, std::move(elem));
  } else
  {
    const auto [min, max] = std::minmax(d1, d2);
    auto elem = std::make_unique<double[]>(max);
    for (auto i = 0; i < min; i++)
    {
      elem[i] = this->arrow[i] - v2.arrow[i];
    }
    if (d1 == min)
      for (auto i = min; i < max; i++)
        elem[i] = -v2.arrow[i];
    else
      for (auto i = min; i < max; i++)
        elem[i] = this->arrow[i];
    return vector(max, std::move(elem));
  }
}

vector subtract(vector& v1, const vector& v2)
{
  return v1.subtract(v2);
}

vector vector::scalar(const double s)
{
  const std::size_t d = this->dimension;
  auto elem = std::make_unique<double[]>(d);
  for (auto i = 0; i < d; i++)
    elem[i] = s * arrow[i];

  return vector(d, std::move(elem));
}

vector vector::unit()
{
  const double norm = this->norm();

  if (compare(norm, 0))
    return 1 * *this; // lazy, probably expensive copy

  const double one_over_norm = 1 / norm;
  const std::size_t d = this->dimension;
  auto elem = std::make_unique<double[]>(d);
  for (auto i = 0; i < d; i++)
    elem[i] = this->arrow[i] * one_over_norm;

  vector u(d, std::move(elem));

  return u;
}

bool vector::equals(const vector& v2) const
{
  const auto d1 = this->dimension;
  const auto d2 = v2.dimension;

  if (d1 != d2)
    return false;
  for (auto i = 0; i < d1; i++)
    if (!compare(this->arrow[i], v2.arrow[i]))
      return false;

  return true;
}

vector operator+(vector& v1, const vector& v2) { return add(v1, v2); }

vector operator-(vector& v1, const vector& v2) { return subtract(v1, v2); }

vector operator*(const double d, vector& v) { return v.scalar(d); }

vector operator*(vector& v, const double d) { return v.scalar(d); }

vector operator/(vector& v, const double d)
{
  if (compare(d, 0.))
    throw vector_exception::div_by_zero();
  return v.scalar(1 / d);
}

bool operator==(const vector& v1, const vector& v2) { return v1.equals(v2); }

constexpr bool compare(const double a, const double b)
{
  constexpr double epsilon = std::numeric_limits<double>::epsilon();

  if constexpr (abs(b - a) < epsilon)
    return true;
  else
    return false;
}

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
  using iterator_category = std::contiguous_iterator_tag;

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
  { return iter_ptr != other.iter_ptr; }

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

  friend auto operator<=>(iterator, iterator) = default;

  /*friend double operator-(iterator first, iterator second);

  friend iterator operator+(iterator, int off);

  friend iterator operator-(iterator, int off);

  friend iterator operator+(int off, iterator);*/

private:
  pointer iter_ptr;
};

vector::iterator vector::begin() { return vector::iterator(arrow.get()); }

vector::iterator vector::end()
{ return vector::iterator(arrow.get() + dimension); }

vector::iterator begin(vector& v) { return v.begin(); }

vector::iterator end(vector& v) { return v.end(); }

struct vector::const_iterator
{
public:

  using value_type = const double;
  using difference_type = std::ptrdiff_t;
  using pointer = const double*;
  using reference = const double&;
  using iterator_category = std::contiguous_iterator_tag;

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

  constexpr difference_type operator-(const_iterator other) const noexcept
  { return iter_ptr - other.iter_ptr; }

  friend auto operator<=>(const_iterator, const_iterator) = default;

private:
  pointer iter_ptr;
};

vector::const_iterator vector::begin() const
{
  return vector::const_iterator(arrow.get());
}

vector::const_iterator vector::end() const
{
  return vector::const_iterator(arrow.get() + dimension);
}

vector::const_iterator vector::cbegin() const
{
  return vector::const_iterator(arrow.get());
}

vector::const_iterator vector::cend() const
{
  return vector::const_iterator(arrow.get() + dimension);
}

/*
 *  END iterator implementation
 */

void vector::print() const
{
  std::cout << "(";
  for (const auto i : *this)
    std::cout << i << ", ";
  std::cout << ")" << std::endl;
}

}
