
import std;
import lam.linearalgebra;

// Stress test compile time by instantiating many constexpr calls
constexpr double run_range()
{
  double sum = 0;
  // Repeat many times to stress constant evaluator
  for (int i = 0; i < 500; ++i)
  {
    lam::linalg::vector<double> v1{1.0, 2.0, 3.0};
    lam::linalg::vector<double> v2{2.0, 3.0, 4.0};
    // Use range algorithms specifically
    sum += lam::linalg::dot_range(v1, v2);
    sum += lam::linalg::distance_range(v1, v2);
    sum += lam::linalg::norm_range(v1);
  }
  return sum;
}

static_assert(run_range() > 0);

int main() { return 0; }
