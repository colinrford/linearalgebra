
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
};

struct vector_info;

class vector {

  private:

    std::size_t dimension;
    std::unique_ptr<double[]> arrow;
    std::optional<std::string> label;
    std::optional<std::unique_ptr<vector_info>> info;

  public:

    vector() = delete; // for now?

    vector(const std::size_t dim);

    vector(const std::size_t dim, std::unique_ptr<double[]> elem);

    vector(const std::vector<double> elem);

    vector(vector&& v) noexcept;

    vector& operator=(vector&& v);

    //const auto operator[](const std::size_t index) const -> decltype(arrow[index]);

    auto operator[](const std::size_t index) -> decltype(arrow[index]);

    [[nodiscard]] constexpr bool empty() const noexcept; // copying stdlib

    constexpr std::size_t get_dimension() const;

    constexpr std::size_t size() const;

    constexpr std::size_t length() const;

    double* get_arrow() const;

    double mag() const;

    double mag2() const;

    double norm() const;

    double norm2() const;

    double dot(const vector& v2) const;

    vector cross(const vector& v2);

    vector add(const vector& v2);

    vector subtract(const vector& v2);

    vector scalar(const double s);

    vector unit();

    bool equals(const vector& v2) const;

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
  double trick;
  double tat;
};

vector add(vector& v1, vector& v2);
vector operator+(vector&, const vector&);
vector operator-(vector&, const vector&);
vector operator*(const double, vector&);
vector operator*(vector&, const double);
vector operator/(vector&, const double);
bool operator==(const vector&, const vector&);
constexpr bool compare(const double, const double);

}
