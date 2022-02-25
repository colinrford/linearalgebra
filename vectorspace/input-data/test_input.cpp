
#include "vs-input.cpp"
#include <ranges>

int main()
{
  auto realm = linalg::data_input::getTerminalInput("directory of data files pls");
  auto possibilities = linalg::data_input::getDataFrom(realm);
  if (!possibilities.has_value())
  {
    std::cout << "no possibilities in this realm" << std::endl;
    return 0;
  }
  auto dataPear = std::move(possibilities.value());
  auto vectorData = std::move(std::get<std::vector<linalg::Vector>>(dataPear));
  auto matrixData = std::move(std::get<std::vector<linalg::Matrix>>(dataPear));

  for (auto& v : vectorData)
    v.print();
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

    std::cout << "has crout determinant " << m.croutLUDet() << std::endl;
    std::cout << "has crout inverse ";
    m_inv.print();
    std::cout << std::endl;
    std::cout << "has doolittle determinant " << m.doolittleLUDet()
              << std::endl;
  }

  return 0;
}
