
#include <vector>
#include <tuple>
#include "vector.h"

// pls compile w at least c++11 

struct MatrixException : public std::exception {
	const char* nonPos () const throw () {
  	return "The number of rows and columns of the matrix must be at least 1";
  }
  const char* addNotDef() const throw () {
  	return "Addition is undefined for these matrices. The numbers of rows and columns must be equal, i.e., n x m + n x m";
  }
  const char* subNotDef() const throw () {
  	return "Subtraction is undefined for these matrices. The numbers of rows and columns must be equal, i.e., n x m - n x m";
  }
  const char* multNotDef() const throw () {
  	return "Multiplication is undefined for these matrices. Matrices must be of the form n x m * m x p";
  }
  const char* notInv() const throw () {
  	return "This matrix is not invertible. Check: Is it square? Is its determinant nonzero?";
  }
};

class Matrix {

	private: 

		int nRows;
		int mColumns;
		std::vector<std::vector<double> > matrix;

	public:

		Matrix(int n, int m);

		Matrix(int n, int m, std::vector<std::vector<double> > mtrx);

		Matrix(Matrix&& mtrx);

		Matrix& operator=(Matrix&& mtrx);

		auto operator[](int index) -> decltype(matrix[index]);

		int getNumRows();

		void setNumRows(int n);

		int getNumColumns();

		void setNumColumns(int m);

		double getEntry(int i, int j);

		std::vector<std::vector<double> > getMatrix();

		void setMatrix(std::vector<std::vector<double> > mtrx);

		Matrix add(Matrix& m_2);

		Matrix subtract(Matrix& m_2);

		Matrix multiply(Matrix& mp);

		Matrix inverse();

		Matrix transpose();

		double determinant();

		double det();

		bool isSquare();

		bool equals(Matrix& m_2);

		void print();

};

Matrix operator+(Matrix&, Matrix&);
Matrix operator-(Matrix&, Matrix&);
Matrix operator*(Matrix&, Matrix&);
Matrix operator*(double, Matrix&);
Matrix operator*(Matrix&, double);
Matrix operator/(Matrix&, double);
bool operator==(Matrix&, Matrix&);
void print(const Matrix&);
