
#include <print>
#include <valarray>
#include "vector_concept.hpp"
#include "../../vectorspace/vector.hpp"

int main()
{
  using namespace linalg;
  static_assert(concepts::vector_c<vector>);
  const vector u{1., 0., 0.};
  const vector v{0., 1., 0.};
  auto lambda = [](concepts::vector_c auto& a, concepts::vector_c auto& b) {
    return a + b;
  };
  const vector w = lambda(u, v);
  w.print();
  linalg::print(w + w);
  std::valarray a{1., 0., 1.};
  std::valarray b{-1., 1., 1.};
  static_assert(concepts::group_element_c_weak<vector>);
  //std::valarray c = lambda(a, b);
  //std::println(c);
}
