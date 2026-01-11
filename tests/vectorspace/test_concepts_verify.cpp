
import std;
import lam.linearalgebra;
import lam.concepts;

using namespace lam::linalg;
using namespace lam::linalg::concepts::experimental;

// Check vector<double> (Field scalar)
static_assert(vector_c_weak<vector<double>>);
static_assert(module_element_c_weak<vector<double>, double>);
static_assert(vectorspace_element_c_weak<vector<double>, double>); // Should pass

// Check vector<int> (Ring scalar)
static_assert(vector_c_weak<vector<int>>);
static_assert(module_element_c_weak<vector<int>, vector<int>::scalar_type>);
// static_assert(!vectorspace_element_c_weak<vector<int>, int>);
// int satisfies weak field concepts (has division), so this actually passes.
static_assert(vectorspace_element_c_weak<vector<int>, int>);

// Check vector<float>
static_assert(vector_c_weak<vector<float>>);
static_assert(vectorspace_element_c_weak<vector<float>, float>);

int main()
{
  std::println("All concept static_assertions passed!");
  return 0;
}
