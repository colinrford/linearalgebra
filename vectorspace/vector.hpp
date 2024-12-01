
#include <cmath>
#include <initializer_list>
#include <iostream>
#include <exception>
#include <limits>
#include <list>
#include <memory>
#include <numeric>
#include <string>
#include <optional>
#include <vector>

namespace linalg
{

// not sure this is the way
struct vector_exception : public std::exception {

  struct non_pos : public std::exception {
    const char* what () const throw () {
      return "The dimension of the vector must be at least 1.\n";
    }
  };
  struct div_by_zero : public std::exception {
    const char* what () const throw () {
      return "Division by zero, ' 0 ', is undefined!\n";
    }
  };
  struct cross_undef : public std::exception {
    const char* what () const throw () {
      return "This cross product is not defined in dimensions other than 3. \
              Wedges later.\n";
    }
  };
  struct out_of_bounds : public std::exception {
    const char* what () const throw () {
      return "Index exceeds dimension of vector.\n";
    }
  };
  struct null_construction : public std::exception {
    const char* what () const throw () {
      return "Construction of a vector with nullptr is not allowed.\n";
    }
  };
  struct plus_equals_unequal_dim : public std::exception {
    const char* what () const throw () {
      return "operator+= is more restricted than others, it requires the \
              vectors have equal dimension.\n";
    }
  };
};

struct vector_info;

class vector {

  private:

    std::size_t dimension;
    std::unique_ptr<double[]> arrow;
    std::optional<std::string_view> label;
    std::optional<std::unique_ptr<vector_info>> info;

    //dangerous
    vector(const std::size_t dim, std::unique_ptr<double[]> elem);

  public:

    using value_type = double;

    vector(const std::size_t);

    vector(const std::initializer_list<value_type>);

    vector(const std::vector<value_type>);

    vector(vector&&) noexcept;

    vector& operator=(vector&&) noexcept;

    //const auto operator[](const std::size_t index) const -> decltype(arrow[index]);

    auto operator[](const std::size_t index) -> decltype(arrow[index]);
    auto operator[](const std::size_t index) const -> decltype(arrow[index]);
    /*
    auto operator[](const int index) -> decltype(arrow[index]);
    auto operator[](const int index) const -> decltype(arrow[index]);
    */
    //[[nodiscard]] constexpr bool empty() const noexcept; // copying stdlib

    constexpr std::size_t get_dimension() const { return dimension; }
    constexpr std::size_t size() const { return dimension; }
    constexpr std::size_t length() const { return dimension; }

    void give_label(std::string);
    void give_label(std::string_view);

    std::string_view get_label() const;
    std::string_view name() const;

    void compute_info();

    value_type* get() const;
    value_type* get_arrow() const;

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

    void print() const;

    ~vector();

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

struct vector_info
{
  double norm;
  double norm_squared;
  constexpr vector_info(double n, double n2) : norm(n), norm_squared(n2)
  { }
};

vector make_zero_vector(const std::size_t);

vector operator+(const vector&, const vector&);
vector operator-(const vector&, const vector&);
vector operator-(const vector&);
vector operator*(const double, const vector&);
vector operator*(const vector&, const double);
vector operator/(const vector&, const double);
bool operator==(const vector&, const vector&);

void print(const vector&);

} // end namespace linalg
