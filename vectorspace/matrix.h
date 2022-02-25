
#include "vector.h"
#include <cassert>
#include <list>
#include <numeric>
#include <utility>

namespace linalg
{

struct LUdcmp;

// also stale
struct MatrixException : public std::exception {
	const char* nonPos () const throw () {
  	return "The number of rows and columns of the matrix must be at least 1";
  }
  const char* addNotDef() const throw () {
  	return "Addition is undefined for these matrices. The numbers of rows \
		 				and columns must be equal, i.e., n x m + n x m";
  }
  const char* subNotDef() const throw () {
  	return "Subtraction is undefined for these matrices. The numbers of \
						rows and columns must be equal, i.e., n x m - n x m";
  }
  const char* multNotDef() const throw () {
  	return "Multiplication is undefined for these matrices. Matrices must be \
						of the form n x m * m x p";
  }
  const char* notInv() const throw () {
  	return "This matrix is not invertible. Check: Is it square? \
						Is its determinant nonzero?";
  }
};


class Matrix {

	private:

		int nRows;
		int mColumns;
		std::optional<std::string> label;
		std::unique_ptr<std::unique_ptr<double[]>[]> matrix;
		std::optional<std::unique_ptr<LUdcmp>> lud;

	public:

		Matrix(int n, int m);

		Matrix(int n, int m, std::unique_ptr<std::unique_ptr<double[]>[]> mtrx);

		Matrix(int n);

		Matrix(int n, std::unique_ptr<std::unique_ptr<double[]>[]> mtrx);

		Matrix(Vector diag);

		Matrix(std::vector<double> diag);

		Matrix(std::vector<std::vector<double> > mtrx);

		Matrix(std::vector<linalg::Vector> rows);

		Matrix(Matrix&& mtrx) noexcept;

		Matrix& operator=(Matrix&& mtrx);

		auto operator[](int index) -> decltype(matrix[index]);

		int getNumRows();

		void setNumRows(int n);

		int getNumColumns();

		void setNumColumns(int m);

		double getEntry(int i, int j);

		//double** getMatrix();

		void setMatrix(std::unique_ptr<std::unique_ptr<double[]>[]> mtrx);

		Matrix add(Matrix& m_2);

		Matrix subtract(Matrix& m_2);

		Matrix multiply(Matrix& mp);

		Matrix inverse();

		Matrix transpose();

		std::optional<linalg::LUdcmp> croutLU();

		Vector croutLUSolveSystem(Vector&);

		Matrix croutLUSolveMatrixSystem(Matrix&);

		double croutLUDet();

		Matrix croutLUInv();

		std::optional<linalg::LUdcmp> doolittleLU();

		double doolittleLUDet();

		double determinant();

		double det();

		bool isSquare();

		bool equals(Matrix& m_2);

		void print();

		void writeTeXto(std::string);

		~Matrix() = default;

};

struct LUdcmp
{
	Matrix lu_decomp;
	Vector permutations;
	int parity;
};

Matrix operator+(Matrix&, Matrix&);
Matrix operator-(Matrix&, Matrix&);
Matrix operator*(Matrix&, Matrix&);
Matrix operator*(double, Matrix&);
Matrix operator*(Matrix&, double);
Vector operator*(Matrix&, Vector);
Matrix operator/(Matrix&, double);
bool operator==(Matrix&, Matrix&);
void print(const Matrix&);

}
