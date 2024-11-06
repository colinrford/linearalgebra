
#include "matrix.cpp"

int main()
{
  auto m = linalg::matrix(4);
  m[1,3] = 2.;
  m[3,1] = 2.;
  m[0,1] = 2.;

  for (auto b = m.column_major_begin(), e = m.column_major_end(); b != e; ++b)
  {
    std::cout << *b << std::endl;
  }
}
