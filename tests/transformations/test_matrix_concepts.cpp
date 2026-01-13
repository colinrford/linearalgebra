/*
 *  test_matrix_concepts.cpp - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  test_matrix_concepts is a c++ module
 */

import std;
import lam.linearalgebra;
import lam.concepts;

using namespace lam::linalg;

// Verify matrix satisfies vectorspace_element_c_weak
static_assert(lam::concepts::experimental::vectorspace_element_c_weak<matrix<double>, double>);
static_assert(lam::concepts::experimental::vectorspace_element_c_weak<matrix<float>, float>);
static_assert(lam::concepts::experimental::vectorspace_element_c_weak<matrix<int>, int>);

int main()
{
  std::println("Matrix satisfies vectorspace_element_c_weak concept.");
  return 0;
}
