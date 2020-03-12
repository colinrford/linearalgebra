#include "matrix.cpp"

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
			throw MatrixException();
	} catch (MatrixException& e) {
		std::cout << "Error! :0( " << e.nonPos() << std::endl;
		return 0;
	}

	std::vector<std::vector<double> > entries(numRows, std::vector<double>(numColumns));
	for (int i = 0; i < numRows; i++) 
	{
		for (int j = 0; j < numColumns; j++)
		{
			std::cout << "Enter value for entry (" << i << ", " << j << ") in matrix 1: ";
			std::cin >> entries[i][j];
		}
	}
	
	std::cout << "\nMatrix m_1 = " << std::endl;
	Matrix m1(numRows, numColumns, std::move(entries));
	m1.print();

	std::vector<std::vector<double> > entries2(numRows, std::vector<double>(numColumns));
	for (int i = 0; i < numRows; i++) 
	{
		for (int j = 0; j < numColumns; j++)
		{
			std::cout << "Enter value for entry (" << i << ", " << j << ") in matrix 2: ";
			std::cin >> entries2[i][j];
		}
	}

	std::cout << "\nMatrix m_2 = " << std::endl;
	Matrix m2(numRows, numColumns, std::move(entries2));
	m2.print();

	std::cout << "\nThe sum of these two matrices is the matrix: " << std::endl;
	Matrix m3 = m1 + m2;
	m3.print();

	std::cout << "\nThe difference of these two matrices is the matrix: " << std::endl;
	Matrix m4 = m1 - m2;
	m4.print();

	bool tru = m1.isSquare();
	if (tru)
		std::cout << "\n~ ~ ~ indeed square ~ ~ ~" << std::endl;
	else
		std::cout << "\nNot At All Square (NAAS) ! ! !" << std::endl;

	std::cout <<"\nThe product of m1 and m2 is the matrix: " << std::endl;
	Matrix m5 = m1 * m2;
	m5.print();

	std::cout << "Hope it worked :0)\nPress return to end program.";
	std::getchar();
	std::getchar();
}