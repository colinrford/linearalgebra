
#include "vector.cpp"
#include <iomanip>
#include <ranges>

template<typename T>
using std_vector = typename std::vector<T>;
using la_vector = typename linalg::vector;

int main()
{
  std::vector wasteful_allocation{1., 2., 3., 4.};
  la_vector another_one(wasteful_allocation);
  auto indexing_set = std::views::iota(1, 100);
  for (auto index : indexing_set)
  {
    auto lp_norm_value = another_one.lpnorm(index);
    std::cout << std::setprecision(22)
              << lp_norm_value << std::endl;
  }
}
