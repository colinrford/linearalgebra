#include "vector.cpp"

int main()
{
  int dim = 3;
  double elem[dim];

  elem[0] = 1;
  elem[1] = -1;
  elem[2] = 2;

  Vector* v1 = new Vector(dim, elem);
  
  v1->print();

  delete v1;
  
  v1 = nullptr;

  return 0;
}


