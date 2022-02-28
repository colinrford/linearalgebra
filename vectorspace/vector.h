
#include <cmath>
#include <iostream>
#include <exception>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
#include <vector>

namespace linalg
{

// stale
struct VectorException : public std::exception {
   const char* nonPos () const throw () {
      return "The dimension of the vector must be at least 1";
   }
   const char* divByZero () const throw () {
      return "Division by zero, ' 0 ', is undefined";
   }
   const char* crossUndef () const throw () {
      return "The cross product is undefined in dimensions other than 3";
   }
};

class Vector {

  private:

    int dimension;
    std::unique_ptr<double[]> arrow;
    std::optional<std::string> label;

  public:

    Vector(int dim);

    Vector(int dim, std::unique_ptr<double[]> elem);

    Vector(std::vector<double> elem);

    Vector(Vector&& v) noexcept;

    Vector& operator=(Vector&& v);

    auto operator[](int index) -> decltype(arrow[index]);

    int getDimension();

    int size();

    int length();

    double* getArrow();

    void setArrow(std::unique_ptr<double[]> elem);

    double norm(Vector& v);

    double dot(Vector& v2);

    Vector cross(Vector& v2);

    Vector add(Vector& v2);

    Vector subtract(Vector& v2);

    Vector scalar(double s);

    Vector unit();

    bool equals(Vector& v2);

    void print() const;

    ~Vector() = default;
};

Vector add(Vector& v1, Vector& v2);
Vector operator+(Vector&, Vector&);
Vector operator-(Vector&, Vector&);
Vector operator*(double, Vector&);
Vector operator*(Vector&, double);
Vector operator/(Vector&, double);
bool operator==(Vector&, Vector&);
bool compare(double, double);

}
