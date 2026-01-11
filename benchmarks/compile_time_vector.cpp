import lam.linearalgebra;
import std;

using namespace lam::linalg;

// Force heavy template instantiation
template<typename T>
void instantiate_stuff()
{
  vector<T> v1(100);
  vector<T> v2(100);
  v1 += v2;
  v1 -= v2;
  v1 *= T{2};
  auto n = v1.norm();
  auto c = v1.cross(v2); // only for T=double/float usually but ok generic
}

void inst()
{
  instantiate_stuff<double>();
  instantiate_stuff<float>();
  instantiate_stuff<int>();
  instantiate_stuff<long>();
}

int main()
{
  inst();
  return 0;
}
