
#include "vector.cpp"

int main()
{
  int dim;
  std::cout << "Enter integer for dimension: ";
  std::cin >> dim;

  try
  {
    if (dim <= 0)
      throw VectorException();
  } catch (VectorException& e) {
      std::cout << "Error! :0( " << e.nonPos() << std::endl;
      return 0;
  }

  unique_ptr<double[]> elem(new double[dim]);
  for (int i = 0; i < dim; i++)
  {
    std::cout << "Enter value for v_1 component " << i << ": ";
    std::cin >> elem[i];
    
  }
  Vector v1(dim, std::move(elem));
  v1.print();

  double norm = v1.norm(v1);
  std::cout << "norm of v_1 = " << norm << "\n";
  
  unique_ptr<double[]> elem2(new double[dim]);
  for (int i = 0; i < dim; i++)
  {
    std::cout << "Enter value for v_2 component " << i << ": ";
    std::cin >> elem2[i];
  }
  Vector v2(dim, std::move(elem2));
  v2.print();

  double norm2 = v2.norm(v2);
  std::cout << "norm of v_2 = " << norm2 << "\n";

  double dp = dot(v1, v2);
  std::cout << "v_1 dot v_2 = " << dp << "\n";

  Vector n(dim);
  try {
    n = cross(v1, v2);
  } catch (VectorException& e) {
    std::cout << "Error! :0( " << e.crossUndef() << std::endl;
    return 0;
  }

  std::cout << "v_1 cross v_2 = ";
  n.print();

  Vector v3 = v1 + v2;
  std::cout << "v_1 + v_2 = ";
  v3.print();
  
  Vector v4 = v1 - v2;
  std::cout << "v_1 - v_2 = ";
  v4.print();

  bool tf = v1 == v2;
  std::cout << "does v_1 = v_2? ";
  if (tf)
    std::cout << "yes\n";
  else
    std::cout << "no\n";

  v1 = 2 * v1;  
  std::cout << "2 * v_1 = "; 
  v1.print();
  v1 = v1 / 2;
  std::cout << "v_1 / 2 = ";
  v1.print();


  std::cout << "Hope it worked :-)\nPress any key to end program.";
  std::getchar();
}
