
module;

import "gcc_headers.hpp";
import la.concepts;
import la.vectors.exceptions;

export module la.vectors;

export namespace la
{
namespace vectors
{
template<typename G>
concept group_element_c_weak = la::concepts::group_element_c_weak<G>;
template<typename K>
concept field_element_c_weak = la::concepts::field_element_c_weak<K>;

template<field_element_c_weak K>
struct vector_info;

template<field_element_c_weak K>
struct vector_info;

template<field_element_c_weak K>
class vector
{
  using size_type = std::size_t;
  using value_type = K;

  size_type dimension;
  std::unique_ptr<value_type[]> arrow;
  std::optional<std::unique_ptr<vector_info<value_type>>> info;

  //dangerous
  vector(const size_type dim, std::unique_ptr<value_type[]> elem);
  value_type* get() const;
  value_type* get_arrow() const;

  public:

  constexpr vector();

  vector(const size_type);

  vector(const std::initializer_list<value_type>);

  vector(const std::vector<value_type>);

  vector(vector&&) noexcept;

  vector& operator=(vector&&) noexcept;

  //const auto operator[](const size_type index) const -> decltype(arrow[index]);

  auto operator[](const size_type index) -> decltype(arrow[index]);
  auto operator[](const size_type index) const -> decltype(arrow[index]);

  [[nodiscard]] constexpr bool empty() const noexcept; // copying stdlib

  constexpr size_type get_dimension() const { return dimension; }
  constexpr size_type size() const { return dimension; }
  constexpr size_type length() const { return dimension; }

  std::string_view get_label() const;
  std::string_view name() const;

  void compute_info();

  value_type dot(const vector&) const;
  value_type better_dot(const vector&) const;
  value_type even_better_dot(const vector&) const;

  value_type mag();
  value_type mag() const;

  value_type mag2();
  value_type mag2() const;

  value_type norm();
  value_type norm() const;

  value_type norm2();
  value_type norm2() const;

  value_type pnorm(const value_type) const;
  value_type lpnorm(const value_type) const;

  vector cross(const vector&) const;

  vector add(const vector&) const;
  vector subtract(const vector&) const;
  vector negate() const;
  vector scalar(const value_type) const;

  vector unit() const;

  constexpr bool equals(const vector&) const;

  vector& operator+=(const vector&);
  vector operator-() const;

  void print() const;

  constexpr ~vector();

  struct iterator;
  iterator begin();
  iterator end();

  struct const_iterator;
  const_iterator begin() const;
  const_iterator end() const;
  const const_iterator cbegin() const;
  const const_iterator cend() const;

  struct unsafe_iterator;
  unsafe_iterator ubegin();
  unsafe_iterator uend();
};  // end class vector

template<field_element_c_weak K>
struct vector_info
{
  K norm;
  K norm_squared;
  constexpr vector_info(K n, K n2) : norm(n), norm_squared(n2)
  { }
};

template<field_element_c_weak K>
vector<K> make_zero_vector(const std::size_t);

template<field_element_c_weak K>
vector<K> operator+(const vector<K>&, const vector<K>&);
template<field_element_c_weak K>
vector<K> operator-(const vector<K>&, const vector<K>&);
template<field_element_c_weak K>
//vector operator-(const vector&);
vector<K> operator*(const K, const vector<K>&);
template<field_element_c_weak K>
vector<K> operator*(const vector<K>&, const K);
template<field_element_c_weak K>
vector<K> operator/(const vector<K>&, const K);
template<field_element_c_weak K>
bool operator==(const vector<K>&, const vector<K>&);

template<field_element_c_weak K>
void print(const vector<K>&);

template<field_element_c_weak K>
constexpr vector<K>::vector() : dimension{0}, arrow{nullptr} {}

template<field_element_c_weak K>
vector<K>::vector(const std::size_t dim)
          : dimension{dim}
{
  if (dim > 0) [[likely]]
  {
    arrow = std::make_unique<K[]>(dim);
    //std::fill(arrow.get(), arrow.get() + dim, 0.);
  } else  throw vector_exception::non_pos();
}

// dangerous, now a private constructor
template<field_element_c_weak K>
vector<K>::vector(const std::size_t dim, std::unique_ptr<K[]> elem)
          : dimension{dim}
{
  if (elem == nullptr) [[unlikely]]
    throw vector_exception::null_construction();
  if (dim == 0) [[unlikely]]
    throw vector_exception::non_pos();
  arrow = std::move(elem);
}

template<field_element_c_weak K>
vector<K>::vector(const std::vector<K> elem)
      : dimension{elem.size()}
{
  if (elem.size() > 0) [[likely]]
  {
    arrow = std::make_unique<K[]>(dimension);
    for (std::size_t i = 0; i < dimension; i++)
      arrow[i] = elem[i];
  } else  throw vector_exception::non_pos();
}

template<field_element_c_weak K>
vector<K>::vector(const std::initializer_list<K> il)
      : dimension{il.size()}
{
  if (il.size() == 0) [[unlikely]]
    throw vector_exception::non_pos();
  arrow = std::make_unique<K[]>(il.size());
  auto arrow_ptr = arrow.get();
  for (auto entry : il)
  {
    *arrow_ptr = entry;
    arrow_ptr++;
  }
}

template<field_element_c_weak K>
vector<K>::vector(vector<K>&& v) noexcept
      : dimension{v.dimension},
        arrow{std::move_if_noexcept(v.arrow)}
{ v.dimension = 0; }

template<field_element_c_weak K>
vector<K>& vector<K>::operator=(vector<K>&& v) noexcept
{
  this->arrow = std::move(v.arrow);
  this->dimension = v.dimension;
  v.dimension = 0;
  return *this;
}

template<field_element_c_weak K>
constexpr bool vector<K>::empty() const noexcept
{
  if (this->dimension == 0)
    return true;
  else
    return false;
}

template<field_element_c_weak K>
constexpr vector<K>::~vector() = default;

} // end namespace vector
template<la::concepts::field_element_c_weak K>
using vector = la::vectors::vector<K>;
} // end namespace la
