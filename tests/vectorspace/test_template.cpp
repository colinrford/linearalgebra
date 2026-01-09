#include <list>
#include <numeric>
#include <sstream>
#include <string>
#include "vectorTemplate.hpp"

int main()
{
  int dim = 3;
  bool awaitingInput = true;
  std::vector<Vector<double>> vectors;
  Vector<double> v(3);
  v[0] = 3.14;
  v[1] = 2.0;
  v[2] = 3.0;
  auto v2 = 2 * v;
  auto v3 = v - v2;
  auto blah = v.dot(v3);
  vectors.push_back(std::move(v));
  vectors.push_back(std::move(v2));
  vectors.push_back(std::move(v3));

  while (awaitingInput)
  {
    //std::cout<< "v dot v3" << v.dot(v3) << std::endl;
    awaitingInput = false;
  }
  std::cout << vectors[2][0] << std::endl;
  std::cout << blah << std::endl;
  auto mid = multiplicative_identity<double>();
  double horses = mid * 3.0;
  std::cout << horses << std::endl;
  auto aid = additive_identity<double>();
  double grapes = aid + 3.0 + aid;
  std::cout << grapes << std::endl;
}
