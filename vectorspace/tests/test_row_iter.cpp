
#include <algorithm>
#include <random>
#include "matrix.cpp"

int main()
{
  auto m = linalg::matrix(4);
  m[1,3] = 2.;
  m[2,1] = 2.;

  auto m2 = linalg::matrix(4);
  m2[1,1] = 3.;
  m2[3,2] = 7.;
  /*std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle(m.row_begin(), m.row_end(), g);*/
  std::sort(m.row_begin(), m.row_end());

  for (auto b = m.row_begin(), e = m.row_end(); b != e; ++b)
  {
    double* bs = b->get();
    for (auto i = 0; i < m.get_num_columns(); i++)
      std::cout << bs[i] << " ";
    std::cout << std::endl;
  }
}
