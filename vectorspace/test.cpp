#include <list>
#include <numeric>
#include <sstream>
#include <string>
#include "vector.cpp"

std::vector<std::string> getVectorInput()
{
    std::string csvInput;
    std::vector<std::string> components;
    std::cout << "Enter vector as csv, no other characters." << std::endl;
    std::getline(std::cin, csvInput);
    std::stringstream s_stream(csvInput); //create string stream from the string
    while(s_stream.good())
    {
        std::string component;
        std::getline(s_stream, component, ','); //get first string delimited by comma
        components.push_back(component);
    }
    return components;
}

linalg::Vector createVectorFrom(std::vector<std::string> components)
{
  std::unique_ptr<double[]> arrow = std::make_unique<double[]>(components.size());
  std::list<int> ell(components.size());
  std::iota(ell.begin(), ell.end(), 0);
  for (auto i : ell)
  {
    arrow[i] = std::stod(components[i]);
  }
  linalg::Vector vec(components.size(), std::move(arrow));
  return std::move(vec);
}

int main()
{
  int dim;
  bool awaitingInput = true;
  std::vector<linalg::Vector> vectors;

  while (awaitingInput)
  {
    std::cout << "Adding vectors to list..." << std::endl;
    auto components = getVectorInput();
    linalg::Vector v = std::move(createVectorFrom(components));
    vectors.push_back(std::move(v));

    std::cout << "Do you want to add more? y/n ";
    std::string response;
    do
    {
      //std::cin.clear(); // repair instream
      // clear the buffer
      //std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      std::getline(std::cin, response);
    }
    while (/*!std::cin.fail() && */response.find('y') == std::string::npos
                            && response.find('n') == std::string::npos);

    if (response.find('y') != std::string::npos)
      awaitingInput = true;
    else if (response.find('n') != std::string::npos)
      awaitingInput = false;

    std::cin.clear();
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }

  for (const auto &v : vectors)
    v.print();

  std::cout << "Hope it worked :-)\nPress return key to end program.";
  std::getchar();
  std::getchar();
}
