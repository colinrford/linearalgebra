
#include <print>
#include "vector_concept.hpp"
#include "../../vectorspace/vector.hpp"

int main()
{
  using namespace linalg;

  const vector u{1., 0., 0.};
  const vector v{0., 1., 0.};
  auto lambda = [] (concepts::group_element_c_weak auto& a,
                    concepts::group_element_c_weak auto& b) {
    return a + b;
  };
  vector w = lambda(u, v);
  w.print();
  linalg::print(w + w);
  //std::println("{}", w);
  static_assert(concepts::group_element_c_weak<vector>);
  static_assert(concepts::group_element_c_weak<double>);
  static_assert(concepts::field_element_c_weak<double>);
  //std::valarray c = lambda(a, b);
  //std::println(c);
}
