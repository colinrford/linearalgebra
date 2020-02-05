
#include "vector.h"
using std::cout;

Vector::Vector(int dim)
{
  _vector = new vect;
  try 
  {
    if (dim <= 0)
      throw VectorException();
  } catch (VectorException& e) 
  {
      std::cout << "Error! :0( " << e.nonPos() << std::endl;
  }
  _vector->dimension = dim;
  _vector->arrow = new double[dim];

  for (int i = 0; i < dim; i++)
  {
    if (i == 0)
      _vector->arrow[i] = 1;
    else
      _vector->arrow[i] = 0;
  }
}

// For now I'll just deal with real numbers and come back later to update for general fields
Vector::Vector(int dim, const double* elem)
{
  _vector = new vect;
  _vector->dimension = dim;
  _vector->arrow = new double[dim];
  
  for (int i = 0; i < dim; i++)
  {
    _vector->arrow[i] = elem[i];
  }
}

double Vector::norm(Vector* v)
{
  return sqrt(dot(v, v));
}

double Vector::dot(Vector* v1, Vector* v2)
{
  if (v1->_vector->dimension != v2->_vector->dimension)
  {
    std::cout << "Undefined.\n";
    return NAN;
  }

  double dot;
  int dim = v1->_vector->dimension;

  for (int i = 0; i < dim; i++)
  {
    dot += v1->_vector->arrow[i] * v2->_vector->arrow[i];
  }
  
  bool boo = compare(dot, 0);
  if (boo)
    dot = 0;

  return dot;
}

double dot(Vector* v1, Vector* v2)
{
  return v1->dot(v1, v2);
}

Vector Vector::cross(Vector* v1, Vector* v2)
{
  if ((v1->_vector->dimension != v2->_vector->dimension) || (v1->_vector->dimension != 3))
  {
    std::cout << "Undefined.\n";
    return NULL;
  }
 
  int dim = v1->_vector->dimension;
  double* elem = new double[v1->_vector->dimension];
  elem[0] = v1->_vector->arrow[1] * v2->_vector->arrow[2] - v1->_vector->arrow[2] * v2->_vector->arrow[1];
  elem[1] = -(v1->_vector->arrow[0] * v2->_vector->arrow[2] - v1->_vector->arrow[2] * v2->_vector->arrow[0]);
  elem[2] = v1->_vector->arrow[0] * v2->_vector->arrow[1] - v1->_vector->arrow[1] * v2->_vector->arrow[0];
  
  for (int i = 0; i < dim; i++)
    if (compare(elem[i], 0))
      elem[i] = 0;

  Vector n = Vector(dim, elem);
  
  return n;
}

Vector cross(Vector* v1, Vector* v2)
{
  return v1->cross(v1, v2);
}

Vector Vector::add(Vector* v1, Vector* v2)
{
  if (v1->_vector->dimension != v2->_vector->dimension)
  {
    std::cout << "Undefined.\n";
    return NAN; 
  }

  int dim = v1->_vector->dimension;
  double* elem = new double[v1->_vector->dimension];
  for (int i = 0; i < dim; i++)
    elem[i] = v1->_vector->arrow[i] + v2->_vector->arrow[i];

  Vector v1v2 = Vector(dim, elem);

  return v1v2;
}

Vector add(Vector* v1, Vector* v2)
{
  return v1->add(v1, v2);
}

Vector Vector::subtract(Vector* v1, Vector* v2)
{
  if (v1->_vector->dimension != v2->_vector->dimension)
  {
    std::cout << "Undefined.\n";
    return NAN;
  }

  double* elem = new double[v1->_vector->dimension];
  for (int i = 0; i < v1->_vector->dimension; i++)
    elem[i] = v1->_vector->arrow[i] - v2->_vector->arrow[i];

  Vector v1mv2 = Vector(v1->_vector->dimension, elem);

  return v1mv2;
}

Vector subtract(Vector* v1, Vector* v2)
{
  return v1->subtract(v1, v2);
}

Vector* Vector::scalar(double s)
{
  for (int i = 0; i < this->_vector->dimension; i++)
  {
    this->_vector->arrow[i] = s * this->_vector->arrow[i];
  }
  
  return this;
}

Vector* Vector::unit()
{
  double norm = this->norm(this);

  for (int i = 0; i < this->_vector->dimension; i++)
  {
    this->_vector->arrow[i] = this->_vector->arrow[i] / norm;
  }
  
  return this;
}

bool Vector::equals(Vector* v1, Vector* v2)
{
  if (v1->_vector->dimension != v2->_vector->dimension)
  {
    return false; 
  } 

  for (int i = 0; i < v1->_vector->dimension; i++)
  {
    if (v1->_vector->arrow[i] != v2->_vector->arrow[i])
    {
      return false;
    }
  }
  
  return true;
}

Vector operator+(Vector v1, Vector v2)
{
  Vector v3 = add(&v1, &v2);
  return v3;
}

Vector operator-(Vector v1, Vector v2)
{
  Vector v3 = subtract(&v1, &v2);
  return v3;
}

Vector operator*(double d, Vector v)
{
  Vector v3 = *v.scalar(d);
  return v3;
}

Vector operator*(Vector v, double d)
{
  Vector u = *v.scalar(d);
  return u;
}

Vector operator/(double d, Vector v)
{
  double d_inv = 1 / d;
  
  std::cout << "ABBE ABBE ABBE before ";

  Vector u = *v.scalar(d_inv);
  return u;
}

Vector operator/(Vector v, double d)
{
  double d_inv = 1 / d;

  std::cout << "ABBE ABBE ABBE after ";
  
  Vector u = *v.scalar(d_inv);
  return u;
}

bool operator==(Vector v1, Vector v2)
{
  return v1.equals(&v1, &v2);
}

void Vector::print()
{
  std::cout << "(";

  for (int i = 0; i < this->_vector->dimension; i++)
  {
    if (i == this->_vector->dimension - 1)
      std::cout << this->_vector->arrow[i];
    else
      std::cout << this->_vector->arrow[i] << ", ";
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
