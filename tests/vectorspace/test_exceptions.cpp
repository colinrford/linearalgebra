#include <cassert>
#include <iostream>

import lam.linearalgebra;

int main()
{
  using namespace lam::linalg;

  try
  {
    matrix<double> m1(2, 3);
    matrix<double> m2(4, 5); // incompatible
    auto m3 = m1 * m2;
    assert(false && "Should have thrown dim_mismatch");
  }
  catch (const matrix_exception::dim_mismatch&)
  {}

  try
  {
    vector<double> v2(4);
    vector<double> v3(4);
    auto c = cross(v2, v3);
    assert(false && "Should have thrown cross_undef");
  }
  catch (const vector_exception::cross_undef&)
  {}

  try
  {
    vector<double> v(3);
    vector<double> z(3);
    auto p = project(v, z); // onto z, which has norm 0
    assert(false && "Should have thrown div_by_zero");
  }
  catch (const vector_exception::div_by_zero&)
  {}

  try
  {
    vector<double> v(3);
    vector<double> z(3);
    auto a = angle(v, z); // z has norm 0
    assert(false && "Should have thrown div_by_zero");
  }
  catch (const vector_exception::div_by_zero&)
  {}

  std::cout << "All exception tests passed." << std::endl;
  return 0;
}
