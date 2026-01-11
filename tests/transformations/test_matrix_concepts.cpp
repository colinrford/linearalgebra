import std;
import lam.linearalgebra;
import lam.concepts;

using namespace lam::linalg;

// Verify matrix satisfies vectorspace_element_c_weak
static_assert(concepts::experimental::vectorspace_element_c_weak<matrix<double>, double>);
static_assert(concepts::experimental::vectorspace_element_c_weak<matrix<float>, float>);
static_assert(concepts::experimental::vectorspace_element_c_weak<matrix<int>, int>);

int main()
{
  std::println("Matrix satisfies vectorspace_element_c_weak concept.");
  return 0;
}
