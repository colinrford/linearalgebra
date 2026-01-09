
#include "matrix.hpp"
#include <ranges>
#include <sstream>

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

auto makeIndexingSet = [](int n) {
	return std::views::iota(0, n);
};

linalg::vector create_vector_from(std::vector<std::string> components)
{
  auto arrow = std::vector<double>(components.size());
  for (std::size_t i : arrow)
    arrow[i] = std::stod(components[i]);
  linalg::vector vec(arrow);
  return vec;
}
/*
void readCSVMatrixInput()
{
	std::string line;
	while(std::getline(stream, line))
	{
		std::istringstream s(line);
		std::string field;
		while (getline(s, field,','))
		{

		}
	}
}*/

linalg::matrix createMatrixFrom(std::vector<linalg::vector> rows)
{
	return std::move(linalg::matrix(std::move(rows)));
}

int main()
{
	int numRows;
	std::cout << "Enter integer for number of rows: ";
	std::cin >> numRows;
	int numColumns;
	std::cout << "Enter integer for number of columns: ";
	std::cin >> numColumns;

	try {
		if ((numRows <= 0) || (numColumns <= 0))
			throw linalg::matrix_exception();
	} catch (linalg::matrix_exception& e) {
		std::cout << "Error! :0( " << e.what() << std::endl;
		return 0;
	}

	auto entries = std::make_unique<std::unique_ptr<double[]>[]>(numRows);
	//std::unique_ptr<std::unique_ptr<double[]>[]> entries(numRows, std::unique_ptr<double[]>(numColumns));
	for (int i = 0; i < numRows; i++)
	{
		for (int j = 0; j < numColumns; j++)
		{
			std::cout << "Enter value for entry (" << i << ", " << j << ") in matrix 1: ";
			std::cin >> entries[i][j];
		}
	}

	std::cout << "\nMatrix m_1 = " << std::endl;
	linalg::matrix m1(numRows, numColumns, std::move(entries));
	m1.print();

	auto entries2 = std::make_unique<std::unique_ptr<double[]>[]>(numRows);//, std::make_unique<double[]>(mColumns));
	//std::unique_ptr<std::unique_ptr<double[]>[]> entries2(numRows, std::unique_ptr<double[]>(numColumns));
	for (int i = 0; i < numRows; i++)
	{
		for (int j = 0; j < numColumns; j++)
		{
			std::cout << "Enter value for entry (" << i << ", " << j << ") in matrix 2: ";
			std::cin >> entries2[i][j];
		}
	}

	std::cout << "\nMatrix m_2 = " << std::endl;
	linalg::matrix m2(numRows, numColumns, std::move(entries2));
	m2.print();

	std::cout << "\nThe sum of these two matrices is the matrix: " << std::endl;
	linalg::matrix m3 = m1 + m2;
	m3.print();

	std::cout << "\nThe difference of these two matrices is the matrix: " << std::endl;
	linalg::matrix m4 = m1 - m2;
	m4.print();

	bool tru = m1.is_square();
	if (tru)
		std::cout << "\n~ ~ ~ indeed square ~ ~ ~" << std::endl;
	else
		std::cout << "\nNot At All Square (NAAS) ! ! !" << std::endl;

	std::cout <<"\nThe product of m1 and m2 is the matrix: " << std::endl;
	linalg::matrix m5 = m1 * m2;
	m5.print();

	std::cout << "Hope it worked :0)\nPress return to end program.";
	std::getchar();
	std::getchar();
}
