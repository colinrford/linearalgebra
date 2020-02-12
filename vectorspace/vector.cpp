
#include "vector.h"
using std::cout;

Vector::Vector(int dim)
{
  try 
  {
    if (dim <= 0)
      throw VectorException();
  } catch (VectorException& e) 
  {
      cout << "Error! :0( " << e.nonPos() << std::endl;
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
Vector::Vector(int dim, unique_ptr<double[]> elem) : arrow(std::move(elem))
{
  try
  {
    if (dim <= 0)
      throw VectorException();
  } catch (VectorException& e)
  {
    cout << "Error! :0( " << e.nonPos() << std::endl;
  }
  
  dimension = dim;
}

Vector::Vector(Vector&& v) : arrow{std::move(v.arrow)}, dimension{v.dimension} {
  v.dimension = 0;
}

Vector& Vector::operator=(Vector&& v) {
     this->setArrow(std::move(v.getArrow()));
     return *this;
}

unique_ptr<double[]>& Vector::getArrow()
{
  return arrow;
}

void Vector::setArrow(unique_ptr<double[]> elem)
{
  arrow = std::move(elem);
}

double Vector::norm(Vector& v)
{
  return sqrt(dot(v, v));
}

double Vector::dot(Vector& v1, Vector& v2)
{
  if (v1.dimension != v2.dimension)
  {
    std::cout << "Undefined.\n";
    //return NAN;
  }

  double dot;
  int dim = v1.dimension;

  for (int i = 0; i < dim; i++)
  {
    dot += v1.arrow[i] * v2.arrow[i];
  }
  
  bool boo = compare(dot, 0);
  if (boo)
    dot = 0;

  return dot;
}

double dot(Vector& v1, Vector& v2)
{
  return v1.dot(v1, v2);
}

Vector Vector::cross(Vector& v1, Vector& v2)
{
  if ((v1.dimension != v2.dimension) || (v1.dimension != 3))
  {
    throw VectorException();
  }
 
  int dim = v1.dimension;
  unique_ptr<double[]> elem(new double[dim]);
  elem[0] = v1.arrow[1] * v2.arrow[2] - v1.arrow[2] * v2.arrow[1];
  elem[1] = -(v1.arrow[0] * v2.arrow[2] - v1.arrow[2] * v2.arrow[0]);
  elem[2] = v1.arrow[0] * v2.arrow[1] - v1.arrow[1] * v2.arrow[0];
  
  for (int i = 0; i < dim; i++)
    if (compare(elem[i], 0))
      elem[i] = 0;

  Vector n(dim, std::move(elem));
  
  return n;
}

Vector cross(Vector& v1, Vector& v2)
{
  return v1.cross(v1, v2);
}

Vector Vector::add(Vector& v1, Vector& v2)
{
  if (v1.dimension != v2.dimension)
  {
    std::cout << "Undefined.\n";
    //return NAN; 
  }

  int dim = v1.dimension;
  unique_ptr<double[]> elem(new double[dim]);
  for (int i = 0; i < dim; i++)
    elem[i] = v1.arrow[i] + v2.arrow[i];

  Vector v1v2(dim, std::move(elem));

  return v1v2;
}

Vector add(Vector& v1, Vector& v2)
{
  return v1.add(v1, v2);
}

Vector Vector::subtract(Vector& v1, Vector& v2)
{
  if (v1.dimension != v2.dimension)
  {
    std::cout << "Undefined.\n";
    //return NAN;
  }

  unique_ptr<double[]> elem(new double[v1.dimension]);
  for (int i = 0; i < v1.dimension; i++)
    elem[i] = v1.arrow[i] - v2.arrow[i];

  Vector v1mv2(v1.dimension, std::move(elem));

  return v1mv2;
}

Vector subtract(Vector& v1, Vector& v2)
{
  return v1.subtract(v1, v2);
}

Vector Vector::scalar(double s)
{
  int dim = this->dimension;
  unique_ptr<double[]> elem = std::move(this->arrow);
  for (int i = 0; i < dim; i++)
  {
    elem[i] = s * elem[i];
  }

  Vector sv(dim, std::move(elem));

  return std::move(sv);
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

bool Vector::equals(Vector& v1, Vector& v2)
{
  if (v1.dimension != v2.dimension)
  {
    return false; 
  } 

  for (int i = 0; i < v1.dimension; i++)
  {
    if (v1.arrow[i] != v2.arrow[i])
    {
      return false;
    }
  }
  
  return true;
}

Vector operator+(Vector& v1, Vector& v2)
{
  Vector v3 = add(v1, v2);
  return v3;
}

Vector operator-(Vector& v1, Vector& v2)
{
  Vector v3 = subtract(v1, v2);
  return v3;
}

Vector operator*(double d, Vector& v)
{
  Vector v3 = v.scalar(d);
  return v3;
}

Vector operator*(Vector& v, double d)
{
  Vector u = v.scalar(d);
  return u;
}

Vector operator/(Vector& v, double d)
{
  double d_inv = 1 / d;
  
  Vector u = v.scalar(d_inv);
  return u;
}

bool operator==(Vector& v1, Vector& v2)
{
  return v1.equals(v1, v2);
}

void Vector::print()
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
  double epsilon = 1E-40;

  if (abs(b - a) < epsilon)
    return true;
  else
    return false;
}
