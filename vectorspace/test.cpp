
//#include <numeric>
#include <sstream>
//#include <string>
//#include <chrono>
#include "vector.cpp"

std::vector<std::string> get_vector_input()
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

linalg::vector create_vector_from(std::vector<std::string> components)
{
  auto arrow = std::make_unique<double[]>(components.size());
  for (auto i = 0; i < components.size(); i++)
    arrow[i] = std::stod(components[i]);
  linalg::vector vec(components.size(), std::move(arrow));
  return vec;
}

int main()
{
  int dim;
  bool awaitingInput = true;
  std::vector<linalg::vector> vectors;

  while (awaitingInput)
  {
    std::cout << "Adding vectors to list..." << std::endl;
    auto components = get_vector_input();
    try {
      linalg::vector v = std::move(create_vector_from(components));
      vectors.push_back(std::move(v));
    } catch (const linalg::vector_exception& e)
    {
      std::cout << e.what();
    }

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


  std::cout << std::endl << std::endl;

  for (auto& v : vectors)
    v.print();

  for (const auto& v : vectors)
  {
    auto t_0 = std::chrono::system_clock::now();
    std::cout << "(";
    for (const auto entry : v)
      std::cout << entry << ", ";
    std::cout<< ")";
    auto t_1 = std::chrono::system_clock::now();
    auto duration = std::chrono
                        ::duration_cast<std::chrono
                        ::nanoseconds>(t_1 - t_0).count();
    std::cout << "range-based for ~" << duration << "ns" << std::endl;
  }

  std::cout << std::endl << std::endl;

  std::cout << "Hope it worked :-)\nPress return key to end program.";
  std::getchar();
  std::getchar();
}
