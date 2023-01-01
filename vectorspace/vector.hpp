
#include <cmath>
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
      return "This cross product is not defined in dimensions other than 3. Wedges later.\n";
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
      return "operator+= is more restricted than others, it requires the vectors have equal dimension.\n";
    }
  };
};

struct vector_info;

class vector {

  private:

    std::size_t dimension;
    std::unique_ptr<double[]> arrow;
    std::optional<std::string> label;
    std::optional<std::unique_ptr<vector_info>> info;

    //dangerous
    vector(const std::size_t dim, std::unique_ptr<double[]> elem);

  public:

    using value_type = double;

    //vector() = delete; // for now?

    vector(const std::size_t dim);

    vector(const std::vector<double> elem);

    vector(vector&& v) noexcept;

    vector& operator=(vector&& v);

    //const auto operator[](const std::size_t index) const -> decltype(arrow[index]);

    auto operator[](const std::size_t index) -> decltype(arrow[index]);

    auto operator[](const std::size_t index) const -> decltype(arrow[index]);

    //[[nodiscard]] constexpr bool empty() const noexcept; // copying stdlib

    constexpr std::size_t get_dimension() const { return dimension; }

    constexpr std::size_t size() const { return dimension; }

    constexpr std::size_t length() const { return dimension; }

    double* get_arrow() const;

    double mag();

    const double mag() const;

    double mag2();

    const double mag2() const;

    double norm();

    const double norm() const;

    double norm2();

    const double norm2() const;

    double dot(const vector& v2) const;

    vector cross(const vector& v2) const;

    vector add(const vector& v2) const;

    vector subtract(const vector& v2) const;

    vector scalar(const double s) const;

    vector unit() const;

    constexpr bool equals(const vector& v2) const;

    vector& operator+=(const vector&);

    void print() const;

    ~vector() = default;

    struct iterator;

    iterator begin();

    iterator end();

    struct const_iterator;

    const_iterator begin() const;

    const_iterator end() const;

    const_iterator cbegin() const;

    const_iterator cend() const;
};

struct vector_info
{
  double norm;
  double norm_squared;
  vector_info(double n, double n2)
  { //gcc figures this out but clang does not
    norm = n;
    norm_squared = n2;
  }
};


vector operator+(const vector&, const vector&);
vector operator-(const vector&, const vector&);
vector operator*(const double, const vector&);
vector operator*(const vector&, const double);
vector operator/(const vector&, const double);
bool operator==(const vector&, const vector&);

constexpr bool compare(const double a, const double b)
{
  constexpr double epsilon = std::numeric_limits<double>::epsilon();

  if (abs(b - a) < epsilon)
    return true;
  else
    return false;
}

}
