
#include <cstdlib>
#include <cmath>
#include <complex>
#include <iostream>

struct vect {
  double* arrow;
  int dimension;                          // dimension of vector; number of elements
};

class Vector {

  vect* _vector;

  public:

    Vector(int dim);                      // create a unit vector of dimension dim

    Vector(int dim, const double* elem);  // create a vector with elements elem and dimension dim

    double norm();                        // norm of this = v; ||v||

    double dot(Vector* v2);                // this = v1 dot with v2; v1 * v2

    Vector* cross(Vector* v2);              // cross this = v1 and v2; v1 x v2

    Vector* add(Vector* v2);                // add this = v1 and v2; v1 + v2

    Vector* subtract(Vector* v2);           // subtract v2 from this = v1; v1 - v2

    Vector* scalar(int s);                 // multiply by scalar; s * this = v1

    Vector* unit();                        // returns the unit vector of given vector this = v; v / ||v||

    bool equals(Vector* v2);               // determines whether two vectors are equal

    //Vector* operator+(const Vector &v);

    //Vector* operator-(const Vector &v);

    //Vector* operator*(const double s);

    //Vector* operator/(const double s);

    //Vector* operator==(const Vector &v);

    void print();                         // print vector contents (x1, x2, ..., xn, ...)

    ~Vector();
};
