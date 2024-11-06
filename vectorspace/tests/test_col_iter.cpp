
#include "matrix.cpp"

int main()
{
  auto m = linalg::matrix(4);
  m[1,3] = 2.;
  m[2,1] = 2.;

  for (auto b = m.row_begin(), e = m.row_end(); b != e; ++b)
  {
    double* bs = b->get();
    for (auto i = 0; i < m.get_num_columns(); i++)
      std::cout << bs[i] << " ";
    std::cout << std::endl;
  }
}
