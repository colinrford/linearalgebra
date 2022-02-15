
#include "vector.h"

namespace linalg
{

Vector::Vector(int dim)
{
  try {
    if (dim <= 0)
      throw VectorException();
  } catch (VectorException& e) {
      std::cout << "Error! :0( " << e.nonPos() << std::endl;
  }
  dimension = dim;
  arrow = std::make_unique<double[]>(dimension);

  for (int i = 0; i < dimension; i++)
  {
    if (i == 0)
      arrow[i] = 1;
    else
      arrow[i] = 0;
  }
}

// For now I'll just deal with real numbers and come back later to update for general fields
Vector::Vector(int dim, std::unique_ptr<double[]> elem) : arrow(std::move(elem))
{
  try {
    if (dim <= 0)
      throw VectorException();
  } catch (VectorException& e) {
    std::cout << "Error! :0( " << e.nonPos() << std::endl;
  }

  dimension = dim;
}

Vector::Vector(std::vector<double> elem)
{
  int dim = elem.size();
  dimension = dim;
  arrow = std::make_unique<double[]>(dimension);
  for (int i = 0; i < dimension; i++)
    arrow[i] = elem[i];
}

Vector::Vector(Vector&& v) noexcept
      : arrow{std::move(v.arrow)}, dimension{v.dimension}
{
  //v.dimension = 0;
}

Vector& Vector::operator=(Vector&& v)
{
     this->setArrow(std::move(v.arrow));
     this->dimension = v.dimension;
     return *this;
}

auto Vector::operator[](int index) -> decltype(arrow[index])
{
  if (index >= this->dimension)
    throw VectorException();

  return arrow[index];
}

int Vector::getDimension()
{
  return dimension;
}

int Vector::size()
{
  return dimension;
}

int Vector::length()
{
  return dimension;
}

// Raw pointers are okay if they are non-owning
double* Vector::getArrow()
{
  return arrow.get();
}

void Vector::setArrow(std::unique_ptr<double[]> elem)
{
  arrow = std::move(elem);
}

double Vector::norm(Vector& v)
{
  return sqrt(v.dot(v));
}

double Vector::dot(Vector& v2)
{
  int d = this->dimension;
  if (d != v2.dimension)
  {
    std::cout << "Undefined.\n";
    //return NAN;
  }

  double dot = 0.0;

  for (int i = 0; i < d; i++)
  {
    dot += this->arrow[i] * v2.arrow[i];
  }

  bool boo = compare(dot, 0);
  if (boo)
    dot = 0;

  return dot;
}

double dot(Vector& v1, Vector& v2)
{
  return v1.dot(v2);
}

Vector Vector::cross(Vector& v2)
{
  int d = this->dimension;
  if ((d!= v2.dimension) || (d != 3))
  {
    throw VectorException();
  }

  auto elem = std::make_unique<double[]>(d);
  elem[0] = this->arrow[1] * v2.arrow[2] - this->arrow[2] * v2.arrow[1];
  elem[1] = -(this->arrow[0] * v2.arrow[2] - this->arrow[2] * v2.arrow[0]);
  elem[2] = this->arrow[0] * v2.arrow[1] - this->arrow[1] * v2.arrow[0];

  for (int i = 0; i < d; i++)
    if (compare(elem[i], 0))
      elem[i] = 0;

  Vector n(d, std::move(elem));

  return std::move(n);
}

Vector cross(Vector& v1, Vector& v2)
{
  return v1.cross(v2);
}

Vector Vector::add(Vector& v2)
{
  int d = this->dimension;
  if (d != v2.dimension)
  {
    std::cout << "Undefined.\n";
    //return NAN;
  }

  auto elem = std::make_unique<double[]>(d);
  for (int i = 0; i < d; i++)
    elem[i] = this->arrow[i] + v2.arrow[i];

  Vector v1v2(d, std::move(elem));

  return v1v2;
}

Vector add(Vector& v1, Vector& v2)
{
  return std::move(v1.add(v2));
}

Vector Vector::subtract(Vector& v2)
{
  int d = this->getDimension();
  if (d != v2.dimension)
  {
    std::cout << "Undefined.\n";
    //return NAN;
  }

  auto elem = std::make_unique<double[]>(d);
  for (int i = 0; i < d; i++)
    elem[i] = this->arrow[i] - v2.arrow[i];

  Vector v1mv2(d, std::move(elem));

  return v1mv2;
}

Vector subtract(Vector& v1, Vector& v2)
{
  return v1.subtract(v2);
}

// The way this method and unit() are written moves the original data in the vector
// Perhaps try something different
Vector Vector::scalar(double s)
{
  int d = this->dimension;
  auto elem = std::move(this->arrow);
  for (int i = 0; i < d; i++)
  {
    elem[i] = s * elem[i];
  }

  Vector sv(d, std::move(elem));

  return sv;
}

Vector Vector::unit()
{
  double norm = this->norm(*this);

  for (int i = 0; i < this->dimension; i++)
  {
    this->arrow[i] = this->arrow[i] / norm;
  }

  Vector u(this->dimension, std::move(this->arrow));

  this->dimension = 0;

  return u;
}

bool Vector::equals(Vector& v2)
{
  int d = this->dimension;

  if (d != v2.dimension)
    return false;

  for (int i = 0; i < d; i++)
    if (this->arrow[i] != v2.arrow[i])
      return false;

  return true;
}

Vector operator+(Vector& v1, Vector& v2)
{
  Vector v3 = add(v1, v2);
  return std::move(v3);
}

Vector operator-(Vector& v1, Vector& v2)
{
  Vector v3 = subtract(v1, v2);
  return std::move(v3);
}

Vector operator*(double d, Vector& v)
{
  Vector v3 = v.scalar(d);
  return std::move(v3);
}

Vector operator*(Vector& v, double d)
{
  Vector u = v.scalar(d);
  return std::move(u);
}

Vector operator/(Vector& v, double d)
{
  double d_inv = 1 / d;

  Vector u = v.scalar(d_inv);
  return std::move(u);
}

bool operator==(Vector& v1, Vector& v2)
{
  return v1.equals(v2);
}

void Vector::print() const
{
  std::cout << "(";

  for (int i = 0; i < this->dimension; i++)
  {
    if (i == this->dimension - 1)
      std::cout << this->arrow[i];
    else
      std::cout << this->arrow[i] << ", ";
  }

  std::cout << ")\n";
}

bool compare(double a, double b)
{
  double epsilon = std::numeric_limits<double>::epsilon();

  if (abs(b - a) < epsilon)
    return true;
  else
    return false;
}

}
