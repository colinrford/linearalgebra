
#include "../vector.cpp"
#include <iomanip>
#include <ranges>
#include <print>

template<typename T = double>
using std_vector = std::vector<T>;
using la_vector = linalg::vector;

int main()
{
  std_vector wasteful_allocation{1., 2., 3., 4.};
  la_vector another_one(wasteful_allocation);
  auto indexing_set = std::views::iota(1, 100);
  for (auto index : indexing_set)
  {
    auto lp_norm_value = another_one.lpnorm(index);
    std::println("{:.20}", lp_norm_value);
  }
}
