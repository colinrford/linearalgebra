

#include "vector.cpp"
using namespace std;

int main()
{
  int dim;
  cout << "Enter integer for dimension";
  cin >> dim;

  double elem[dim];
  for (int i = 0; i < dim; i++)
  {
    cout << "Enter double value for vector element " << i;
    cin >> elem[i];
  }
  
  Vector v1 = Vector(dim, elem);

  v1.print();

  cout << "Hope it worked :-)";
}
