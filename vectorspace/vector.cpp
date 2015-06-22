
#include "vector.h"

Vector::Vector(int dim)
{
  _vector* = new vect;
  _vector.dimension = dim;
  _vector.arrow = new double[dim];

  for (int i = 0; i < dim; i++)
  {
    if (i == 0)
      _vector.arrow[i] = 1;
    else
      _vector.arrow[i] = 0;
  }
}

// For now I'll just deal with real numbers and come back later to update for general fields
Vector::Vector(int dim, const double* elem)
{
  _vector = new vect;
  _vector.dimension = dim;
  _vector.arrow = new double[dim];
  
  for (int i = 0; i < dim; i++)
  {
    _vector.arrow[i] = elem[i];
  }
}

double Vector::norm()
{
  return sqrt(this.dot(this));
}

double Vector::dot(Vector v2)
{
  if (this._vector.dimension != v2._vector.dimension)
  {
    printf("Undefined.\n");
    return 0;
  }

  double dot;
  int dim = this._vector.dim;

  for (int i = 0; i < dim; i++)
  {
    dot += this._vector.arrow[i] * v2._vector.arrow[i];
  }

  return dot;
}

Vector Vector::cross(Vector v2)
{
  if (this._vector.dimension != v2._vector.dimension || this._vector.dimension != 3)
  {
    printf("Undefined.\n");
    return NULL;
  }
 
  int dim = this._vector.dimension;
  double* elem = new double[this._vector.dimension];
  elem[0] = this._vector.arrow[1] * v2._vector.arrow[2] - this._vector.arrow[2] * v2._vector.arrow[1];
  elem[1] = -(this._vector.arrow[0] * v2._vector.arrow[2] - this._vector.arrow[2] * v2._vector.arrow[0]);
  elem[2] = this._vector.arrow[0] * v2._vector.arrow[1] - this._vector.arrow[1] * v2._vector.arrow[0];
  
  Vector n = new Vector(dim, elem);
  
  return n;
}

Vector Vector::add(Vector v2)
{
  if (this._vector.dimension != v2._vector.dimension)
  {
    printf("Undefined.\n");
  }

}

Vector Vector::subtract(Vector v2)
{
  if (this._vector.dimension != v2._vector.dimension)
  {
    printf("Undefined.\n");
  } 

  return ;
}

Vector Vector::scalar(int s)
{
  for (int i = 0; i < this._vector.dimension; i++)
  {
    this._vector.arrow[i] = s * this._vector.arrow[i];
  }
  
  return this;
}

Vector Vector::unit()
{
  double norm = this.norm();

  for (int i = 0; i < this._vector.dimension; i++)
  {
    this._vector.arrow[i] = this._vector.arrow[i] / norm;
  }
  
  return this;
}

bool Vector::equals(Vector v2)
{
  if (this._vector.dimension != v2._vector.dimension)
  {
    return false; 
  } 

  for (int i = 0; i < this._vector.dimension; i++)
  {
    if (this._vector.arrow[i] != v2._vector.arrow[i])
    {
      return false;
    }
  }
  
  return true;
}

void Vector::print()
{
  cout << "[";

  for (int i = 0; i < this._vector.dimension; i++)
  {
    cout << this._vector.arrow[i] << ", ";
  }

  cout << "]";
}
