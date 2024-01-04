
#include "vs_input.cpp"
//#include <ranges>

int main()
{
  auto realm = "."; // checks current folder
                    //linalg::data_input
                    //      ::get_terminal_input("directory of data files pls");
  auto possibilities = linalg::data_input::get_data_from(realm);
  if (!possibilities.has_value())
  {
    std::cout << "no possibilities in this realm" << std::endl;
    return 0;
  }
  auto dataPear = std::move(possibilities.value());
  auto vectorData = std::move(std::get<std::vector<linalg::vector>>(dataPear));
  auto matrixData = std::move(std::get<std::vector<linalg::matrix>>(dataPear));

  for (auto& v : vectorData)
    v.print();
  std::cout << std::endl;
  for (auto& m : matrixData)
  {
    m.print();
    std::cout << std::endl;
    auto m_squared = m * m;
    m_squared.print();
    std::cout << std::endl;
    auto twice_m = m + m;
    twice_m.print();
    std::cout << std::endl;
    auto m_inv = m.croutLUInv();
    auto m_times_m_inv = m * m_inv;
    auto m_inv_times_m = m_inv * m;
    bool equ = m_times_m_inv == m_inv_times_m;
    linalg::matrix id_mtrx(m.get_num_rows());
    bool equ_ide = m_times_m_inv == id_mtrx;
    if (equ && equ_ide)
      std::cout << "inv seems ok" << std::endl;
    else
      std::cout << "something wrogn w inverse" << std::endl;

    std::cout << "has crout determinant " << m.croutLUDet() << std::endl;
    std::cout << "has crout inverse ";
    m_inv.print();
    std::cout << std::endl;
    std::cout << "has doolittle determinant " << m.doolittleLUDet()
              << std::endl;
  }

  return 0;
}
