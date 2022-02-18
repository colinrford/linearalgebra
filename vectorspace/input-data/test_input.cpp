
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
    m.print();

  return 0;
}
