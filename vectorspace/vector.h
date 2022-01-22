
#include <cmath>
#include <iostream>
#include <exception>
#include <limits>
#include <memory>
#include <vector>

namespace linalg
{

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

  public:

    Vector(int dim);                      // create a unit vector of dimension dim

    Vector(int dim, std::unique_ptr<double[]> elem);  // create a vector with elements elem and dimension dim

    Vector(Vector&& v);

    Vector& operator=(Vector&& v);

    auto operator[](int index) -> decltype(arrow[index]);

    int getDimension();

    int size();

    int length();

    double* getArrow();       // accessor for arrow

    void setArrow(std::unique_ptr<double[]> elem);

    double norm(Vector& v);                        // norm of this = v; ||v||

    double dot(Vector& v2);                // this = v1 dot with v2; v1 * v2

    Vector cross(Vector& v2);              // cross this = v1 and v2; v1 x v2

    Vector add(Vector& v2);                // add this = v1 and v2; v1 + v2

    Vector subtract(Vector& v2);           // subtract v2 from this = v1; v1 - v2

    Vector scalar(double s);                 // multiply by scalar; s * this = v1

    Vector unit();                        // returns the unit vector of given vector this = v; v / ||v||

    bool equals(Vector& v2);               // determines whether two vectors are equal

    void print() const;                         // print vector contents (x1, x2, ..., xn, ...)

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
