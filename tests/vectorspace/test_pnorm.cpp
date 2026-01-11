import std;
import lam.linearalgebra;

int main()
{
  std::vector wasteful_allocation{1., 2., 3., 4.};
  lam::linalg::vector another_one(wasteful_allocation);

  for (int index = 1; index < 10; ++index)
  {
    auto lp_norm_value = another_one.pnorm((double)index);
    std::println("{:.20}", lp_norm_value);
  }
}
