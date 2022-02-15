
#include "matrix.cpp"
#include <sstream>

linalg::Vector createVectorFrom(std::vector<std::string> components)
{
  auto arrow = std::make_unique<double[]>(components.size());
	auto ell = linalg::makeIndexingSet(components.size());

  for (auto i : ell)
    arrow[i] = std::stod(components[i]);
  linalg::Vector vec(components.size(), std::move(arrow));
  return std::move(vec);
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

linalg::Matrix createMatrixFrom(std::vector<linalg::Vector> rows)
{
	return std::move(linalg::Matrix(std::move(rows)));
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
			throw linalg::MatrixException();
	} catch (linalg::MatrixException& e) {
		std::cout << "Error! :0( " << e.nonPos() << std::endl;
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
	linalg::Matrix m1(numRows, numColumns, std::move(entries));
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
	linalg::Matrix m2(numRows, numColumns, std::move(entries2));
	m2.print();

	std::cout << "\nThe sum of these two matrices is the matrix: " << std::endl;
	linalg::Matrix m3 = m1 + m2;
	m3.print();

	std::cout << "\nThe difference of these two matrices is the matrix: " << std::endl;
	linalg::Matrix m4 = m1 - m2;
	m4.print();

	bool tru = m1.isSquare();
	if (tru)
		std::cout << "\n~ ~ ~ indeed square ~ ~ ~" << std::endl;
	else
		std::cout << "\nNot At All Square (NAAS) ! ! !" << std::endl;

	std::cout <<"\nThe product of m1 and m2 is the matrix: " << std::endl;
	linalg::Matrix m5 = m1 * m2;
	m5.print();

	std::cout << "Hope it worked :0)\nPress return to end program.";
	std::getchar();
	std::getchar();
}
