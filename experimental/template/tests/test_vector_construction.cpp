
import la.concepts;
import la.vectors;
import <print>;

int main()
{
  static_assert(la::concepts::field_element_c_weak<double>);
  constexpr la::vector<double> uve = la::vector<double>();
  if (uve.empty())
    std::println("constexpr trivial construction successful.");
  else
    std::println("constexpr trivial construction failed.");
}
