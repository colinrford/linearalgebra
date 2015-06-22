
#include "vector.cpp"

int main()
{
  int dim;
  std::cout << "Enter integer for dimension: ";
  std::cin >> dim;

  double elem[dim];
  for (int i = 0; i < dim; i++)
  {
    std::cout << "Enter double value for vector element " << i << ": ";
    std::cin >> elem[i];
  }
  
  Vector v1 = Vector(dim, elem);

  v1.print();

  std::cout << "Hope it worked :-)\n";
}
