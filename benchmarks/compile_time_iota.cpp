/*
 *  compile_time_iota.cpp - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  compile_time_iota is a c++ module
 */

import std;
import lam.linearalgebra;

// Define a local helper to test the concept
template<typename V1, typename V2>
constexpr auto dot_iota(const V1& v1, const V2& v2)
{
  using T = typename V1::scalar_type;
  auto n = std::min(v1.size(), v2.size());
  // Use iota view to generate indices
  auto indices = std::views::iota(std::size_t{0}, n);
  // fold_left over indices, accessing elements by index
  return std::ranges::fold_left(indices, T{0}, [&](T sum, auto i) { return sum + v1[i] * v2[i]; });
}

constexpr double run_iota_test()
{
  double sum = 0;
  for (int i = 0; i < 50; ++i)
  {
    // Use vector (dynamic allocation) to test the problematic case
    lam::linalg::vector<double> v1{1.0, 2.0, 3.0};
    lam::linalg::vector<double> v2{2.0, 3.0, 4.0};
    sum += dot_iota(v1, v2);
  }
  return sum;
}

static_assert(run_iota_test() > 0);

int main() { return 0; }
