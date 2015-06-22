
#include "vector.cpp"

int main()
{
  int dim;
  std::cout << "Enter integer for dimension: ";
  std::cin >> dim;

  double elem[dim];
  for (int i = 0; i < dim; i++)
  {
    std::cout << "Enter double value for v1 element " << i << ": ";
    std::cin >> elem[i];
  }
  Vector v1 = Vector(dim, elem);
  v1.print();

  double norm = v1.norm();
  std::cout << "v1 norm = " << norm << "\n";
  
  double elem2[dim];
  for (int i = 0; i < dim; i++)
  {
    std::cout << "Enter double value for v2 element " << i << ": ";
    std::cin >> elem2[i];
  }
  Vector v2 = Vector(dim, elem2);
  v2.print();

  double norm2 = v2.norm();
  std::cout << "v2 norm = " << norm2 << "\n";

  double dp = v1.dot(&v2);
  std::cout << "v1 dot v2 = " << dp << "\n";

  Vector n = v1.cross(&v2);
  std::cout << "v1 cross v2 = ";   
  n.print();

  Vector v1v2 = v1.add(&v2);   
  std::cout << "v1 + v2 = ";
  v1v2.print();

  


  v1 = *v1.scalar(2);  
  std::cout << "2 * v1 = "; 
  v1.print();



  std::cout << "Hope it worked :-)\n";

}
