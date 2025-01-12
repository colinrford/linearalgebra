
#include <print>
#include <string>
#include <vector>

import la.vectors;

int main()
{
  constexpr la::vector uve;
  if constexpr (uve.empty())
    std::println("constexpr trivial construction successful.");
  else
    std::println("constexpr trivial construction failed.");
  {
    try {
      la::vector throwaway{};
      std::vector<double> empty_vector;
      std::println("\tgoing to construct la::vector from empty std::vector...");
      la::vector another_one(empty_vector);
    } catch (la::vector_exception::non_pos& e) {
      const char* error_string = e.what();
      std::print("\tla::vector_exception::non_pos: {}", error_string);
    }
  }
  {
    la::vector u(1);
    if (!u.empty())
      std::println("std::size_t constructor succeeded.");
  }
  {
    std::vector sv{1.,2.,3.};
    la::vector vs(sv);
    if (!vs.empty())
      std::println("std::vector constructor succeeded.");
  }
  {
    la::vector w{1.,2.,3.};
    if (!w.empty())
      std::println("initializer_list constructor succeeded.");
  }
  la::vector initially_empty_vector;
  if (initially_empty_vector.empty())
    std::println("\tinitially_empty_vector.size() = 0...");
  {
    la::vector scoped_vector{1.,2.,3.,4.};
    initially_empty_vector = std::move(scoped_vector);
    std::println("\tmoved");
  }
  if (!initially_empty_vector.empty())
    std::println("move assignment succeeded.");
  else
    std::println("move assignment failed.");
}
