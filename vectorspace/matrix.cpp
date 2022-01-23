
#include "matrix.h"
#include <cassert>
#include <list>
#include <numeric>

namespace linalg
{

auto makeIndexingSet = [](int n) -> std::list<int> {
	std::list<int> ell(n);
	std::iota(ell.begin(), ell.end(), 0);
	return ell;
};

auto makeMatrix_unique = [](int numRows, int numColumns) {
	auto mtrx = std::make_unique<std::unique_ptr<double[]>[]>(numRows);
	auto rose = makeIndexingSet(numRows);
	for (auto rho : rose)
		mtrx[rho] = std::make_unique<double[]>(numColumns);

	return std::move(mtrx);
};

// The basic constructor just creates an n x m identity matrix
Matrix::Matrix(int n, int m)
{
	try {
		if (n <= 0 || m <= 0)
			throw MatrixException();
	} catch (MatrixException& e) {
			std::cout << "Error :0( " << e.nonPos() << std::endl;
	}

	nRows = n;
	mColumns = m;
	auto rose = makeIndexingSet(nRows);
	auto calls = makeIndexingSet(mColumns);

	matrix = makeMatrix_unique(n, m);

	for (auto rho : rose)
		for (auto xi : calls)
			if (rho == xi)
				matrix[rho][xi] = 1;
			else
				matrix[rho][xi] = 0;
}

Matrix::Matrix(int n, int m,
								std::unique_ptr<std::unique_ptr<double[]>[]> mtrx)
								: matrix{std::move(mtrx)}
{
	try {
		if (n <= 0 || m <= 0)
			throw MatrixException();
	} catch (MatrixException& e) {
			std::cout << "Error :0( " << e.nonPos() << std::endl;
	}

	nRows = n;
	mColumns = m;
}

Matrix::Matrix(int n)
{
	try {
		if (n <= 0)
			throw MatrixException();
	} catch (MatrixException& e) {
			std::cout << "Error :0( " << e.nonPos() << std::endl;
	}

	nRows = mColumns = n;
	auto square = makeIndexingSet(nRows);

	matrix = makeMatrix_unique(n, n);

	for (auto rho : square)
		for (auto xi : square)
			if (rho == xi)
				matrix[rho][xi] = 1;
			else
				matrix[rho][xi] = 0;
}

Matrix::Matrix(int n, std::unique_ptr<std::unique_ptr<double[]>[]> mtrx)
								: matrix{std::move(mtrx)}
{
	try {
		if (n <= 0)
			throw MatrixException();
	} catch (MatrixException& e) {
			std::cout << "Error :0( " << e.nonPos() << std::endl;
	}

	nRows = mColumns = n;
}

Matrix::Matrix(std::vector<std::vector<double> > mtrx)
{
	auto rose = makeIndexingSet(mtrx.size());
	auto calls = makeIndexingSet(mtrx.front().size());

	matrix = makeMatrix_unique(rose.size(), calls.size());

	for (auto rho : rose)
		for (auto xi : calls)
			matrix[rho][xi] = mtrx[rho][xi];
}

Matrix::Matrix(std::vector<Vector> rows)
{
	int numColumns = rows.front().getDimension();
	auto rose = makeIndexingSet(rows.size());
	auto calls = makeIndexingSet(numColumns);

	matrix = makeMatrix_unique(rows.size(), calls.size());

	for (auto rho : rose)
		for (auto xi : calls)
			matrix[rho][xi] = rows[rho][xi];
}

Matrix::Matrix(std::vector<double> diag)
{
	int n = diag.size();
	auto square = makeIndexingSet(n);

	matrix = makeMatrix_unique(n, n);

	for (auto rho : square)
		for (auto xi : square)
			if (rho != xi)
				matrix[rho][xi] = 0.0;
			else
				matrix[rho][xi] = diag[rho];
}

Matrix::Matrix(Vector diag)
{
	int n = diag.getDimension();
	auto square = makeIndexingSet(n);

	matrix = makeMatrix_unique(n, n);

	for (auto rho : square)
		for (auto xi : square)
			if (rho != xi)
				matrix[rho][xi] = 0.0;
			else
				matrix[rho][xi] = diag[rho];
}

Matrix::Matrix(Matrix&& mtrx) : matrix{std::move(mtrx.matrix)},
																nRows{mtrx.nRows},
																mColumns{mtrx.mColumns}
{
  mtrx.nRows = mtrx.mColumns = 0;
}


Matrix& Matrix::operator=(Matrix&& mtrx)
{
     this->setMatrix(std::move(mtrx.matrix));
     return *this;
}


auto Matrix::operator[](int index) -> decltype(matrix[index])
{
	int numRows = this->getNumRows();
	if (index >= numRows)
		throw MatrixException();

	return this->matrix[index];
}

/*double** Matrix::getMatrix()
{
	//auto m = matrix.get();

	return this->matrix;
}*/

void Matrix::setMatrix(std::unique_ptr<std::unique_ptr<double[]>[]> mtrx)
{
	matrix = std::move(mtrx);
}

int Matrix::getNumRows()
{
	return this->nRows;
}

int Matrix::getNumColumns()
{
	return this->mColumns;
}

void Matrix::setNumRows(int rows)
{
	this->nRows = rows;
}

void Matrix::setNumColumns(int columns)
{
	this->mColumns = columns;
}

double Matrix::getEntry(int i, int j)
{
	int numRows = this->getNumRows();
	int numColumns = this->getNumColumns();

	if ((i > numRows) || (j > numColumns))
		throw MatrixException();

	return this->matrix[i][j];
}

bool canAdd(Matrix& m_1, Matrix& m_2)
{
	int rowsM_1 = m_1.getNumRows();
	int rowsM_2 = m_2.getNumRows();

	int columnsM_1 = m_1.getNumColumns();
	int columnsM_2 = m_2.getNumColumns();

	if ((rowsM_1 != rowsM_2) || (columnsM_1 != columnsM_2))
		return false;
	else
		return true;
}

Matrix Matrix::add(Matrix& m_2)
{
	if (!canAdd(*this, m_2))
		throw MatrixException();

	int numRows = this->getNumRows();
	int numColumns = this->getNumColumns();
	auto rose = makeIndexingSet(numRows);
	auto calls = makeIndexingSet(numColumns);

	auto entries = makeMatrix_unique(numRows, numColumns);

	for (auto rho : rose)
		for (auto xi : calls)
			entries[rho][xi] = this->matrix[rho][xi] + m_2.matrix[rho][xi];

	Matrix m1pm2(numRows, numColumns, std::move(entries));

	return std::move(m1pm2);
}

Matrix add(Matrix& m_1, Matrix& m_2)
{
	return std::move(m_1.add(m_2));
}

bool canSubtract(Matrix& m_1, Matrix& m_2)
{
	return canAdd(m_1, m_2);
}

Matrix Matrix::subtract(Matrix& m_2)
{
	if (!canSubtract(*this, m_2))
		throw MatrixException();

	int numRows = this->getNumRows();
	int numColumns = this->getNumColumns();
	auto rose = makeIndexingSet(numRows);
	auto calls = makeIndexingSet(numColumns);

	auto entries = makeMatrix_unique(numRows, numColumns);

	for (auto rho : rose)
		for (auto xi : calls)
			entries[rho][xi] = this->matrix[rho][xi] - m_2.matrix[rho][xi];

	Matrix m1pm2(numRows, numColumns, std::move(entries));

	return std::move(m1pm2);
}

Matrix subtract(Matrix& m_1, Matrix& m_2)
{
	return std::move(m_1.subtract(m_2));
}

bool canMultiply(Matrix& m_1, Matrix& m_2)
{
	int columnsM_1 = m_1.getNumColumns();
	int rowsM_2 = m_2.getNumRows();

	if (columnsM_1 == rowsM_2)
		return true;
	else
		return false;
}

Matrix Matrix::multiply(Matrix& mp)
{
	if (!canMultiply(*this, mp))
		throw MatrixException();

	int numRowsNM = this->getNumRows();
	int numColumnsNM = this->getNumColumns();
	int numRowsMP = mp.getNumRows();
	assert(numColumnsNM == numRowsMP);
	int numColumnsMP = mp.getNumColumns();

	auto rowsNP = makeIndexingSet(numRowsNM);
	auto columnsNP = makeIndexingSet(numColumnsMP);
	auto columnsNM = makeIndexingSet(numColumnsNM);

	auto entries = makeMatrix_unique(rowsNP.size(), columnsNP.size());

	for (auto i : rowsNP)
		for (auto j : columnsNP)
			for (auto k : columnsNM)
				entries[i][j] += this->matrix[i][k] * mp.matrix[k][j];

	Matrix np(rowsNP.size(), columnsNP.size(), std::move(entries));

	return std::move(np);
}

Matrix multiply(Matrix& nm, Matrix& mp)
{
	return std::move(nm.multiply(mp));
}

bool isInvertible(Matrix& nm)
{
	if (!nm.isSquare())
	{
		throw MatrixException();
		return false;
	}
	if (compare(0, nm.det()))
		return false;

	return true;
}

// TODO: implement
Matrix Matrix::inverse()
{
	if (!isInvertible(*this))
	{
		throw MatrixException();
		//return false;
	}

	int numRows = this->getNumRows();
	int numColumns = this->getNumColumns();

	auto entries = makeMatrix_unique(numRows, numColumns);

	Matrix m(numRows, numColumns, std::move(entries));

	return std::move(m);
}

void gaussje()
{

}

void gaussjInv()
{

}

// TODO: implement
double Matrix::determinant()
{
	return 12345;
}

double Matrix::det()
{
	return this->determinant();
}

// from Numerical Recipes, 3rd Ed – A is square
// std::optional w Matrix, optionally singular?
// Crout uses unit diagonals for the upper triangle
std::optional<std::pair<Matrix, int>> Matrix::croutLU()
{
	int imax = 0, parity = 1;
	int n = this->getNumRows();
	double big, temp;
	auto square = makeIndexingSet(n);
	Vector implicitScalingPerRow(n);
	auto lu = makeMatrix_unique(n, n);

	for (auto i : square)
	{
		big = 0.0;
		for (auto j : square)
		{
			lu[i][j] = this->matrix[i][j];
			if ((temp = abs(lu[i][j])) > big)
				big = temp;
		}
		if (compare(big, 0.0))
		{
			//throw("singular maatrix in LUdcmp");
			return std::nullopt;
		}
		implicitScalingPerRow[i] = 1.0 / big;
	}
	for (auto k : square)
	{
		big = 0.0;
		for (int i = k; i < n; i++)
		{
			temp = implicitScalingPerRow[i] * abs(lu[i][k]);
			if (temp > big)
			{
				big = temp;
				imax = i;
			}
		}
		if (k != imax)
		{
			for (auto j : square)
			{
				temp = lu[imax][j];
				lu[imax][j] = lu[k][j];
				lu[k][j] = temp;
			}
			parity = -parity;
			implicitScalingPerRow[imax] = implicitScalingPerRow[k];
		}
		if (compare(lu[k][k], 0.0)) // if 0, matrix is singular
			lu[k][k] = std::numeric_limits<double>::epsilon();
		for (int i = k + 1; i < n; i++)
		{
			temp = lu[i][k] /= lu[k][k];
			for (int j = k + 1; j < n; j++)
				lu[i][j] -= temp * lu[k][j];
		}
	}
	Matrix lu_decomposition = Matrix(n, std::move(lu));
	return std::make_pair(std::move(lu_decomposition), parity);
}

double Matrix::croutLUDet()
{
	std::pair<Matrix, int> pairity = this->croutLU().value();
	if (pairity.second == 1 || pairity.second == -1)
	{
		double dtrmnnt = pairity.second;
		auto diag = makeIndexingSet(this->nRows);

		for (auto indx : diag)
			dtrmnnt *= pairity.first[indx][indx];

		return dtrmnnt;
	}
	else
		return 0.0;
}

void Matrix::croutLUSolveSystem()
{
	int i, j, m = b.ncols();
	int n = this->getNumRows();
	if (b.nrows() != n || x.nrows() != n || b.ncols() != x.ncols())
		throw("LUdcmp::solve bad sizes");
	Vector implicitScalingPerRow(n);
	for (j = 0; j < m; j++)
	{
		for (i = 0; i < n; i++)
			xx[i] = b[i][j];
		solve(xx, xx);
		for (i = 0; i < n; i++)
			x[i][j] = xx[i];
	}
}

Matrix Matrix::croutLUInv()
{
	std::pair<Matrix, int> pairity = this->croutLU().value();
	int n = pairity.first.getNumRows();
	auto square = makeIndexingSet(n);
	auto ainv = makeMatrix_unique(n, n);

  for (int i : square)
	{
		for (int j : square)
			ainv[i][j] = 0.0;
    ainv[i][i] = 1.0;
  }
  croutLUSolveSystem(ainv, ainv);
}

// adapted from:	sci.utah.edu/~wallstedt
// main diagonal of L is composed with 1s
Matrix Matrix::doolittleLU()
{
	int n = this->getNumRows();
	auto square = makeIndexingSet(n);
	auto lu = makeMatrix_unique(n, n);
	for (int k : square)
	{
		auto subsq = makeIndexingSet(k);
		auto subsqc = makeIndexingSet(n - k);
		auto subsqcpo = makeIndexingSet(n - k + 1);
    for (int j : subsqc)//(int j = k; j < n; ++j)
		{
      double sum = 0.0;
      for (int p : subsq)
				sum += lu[k][p] * lu[p][j];
      lu[k][j] = (this->matrix[k][j] - sum);
    }
    for (int i : subsqcpo)//(int i = k + 1; i < n; ++i)
		{
    	double sum = 0.0;
    	for (int p : subsq)
				sum += lu[i][p] * lu[p][k];
    	lu[i][k] = (this->matrix[i][k] - sum) / lu[k][k];
    }
	}
	Matrix lu_decomposition(n, std::move(lu));
	return lu_decomposition;
}
//incomplete possibly may need permutation info too?
double Matrix::doolittleLUDet()
{
	Matrix lu_decomposition = this->doolittleLU();

	double dtrmnnt = 0.0;
	auto diag = makeIndexingSet(this->nRows);

	for (auto indx : diag)
		dtrmnnt *= lu_decomposition[indx][indx];

	return dtrmnnt;
}

bool Matrix::isSquare()
{
	int numRows = this->getNumRows();
	int numColumns = this->getNumColumns();

	if (numRows == numColumns)
		return true;
	else
		return false;
}

bool Matrix::equals(Matrix& m_2)
{
	int numRowsM_1 = this->getNumRows();
	int numColumnsM_1 = this->getNumColumns();
	int numRowsM_2 = m_2.getNumRows();
	int numColumnsM_2 = m_2.getNumColumns();

  if ((numRowsM_1 != numRowsM_2) || (numColumnsM_1 != numColumnsM_2))
  {
    return false;
  }

	auto rose = makeIndexingSet(numRowsM_1);
	auto calls = makeIndexingSet(numColumnsM_1);

  for (auto rho : rose)
  	for (auto xi : calls)
    	if (!compare(this->matrix[rho][xi], m_2.matrix[rho][xi]))
     		return false;

  return true;
}

Matrix operator+(Matrix& m_1, Matrix& m_2)
{
  return std::move(add(m_1, m_2));
}

Matrix operator-(Matrix& m_1, Matrix& m_2)
{
  return std::move(subtract(m_1, m_2));
}

Matrix operator*(Matrix& m_1, Matrix& m_2)
{
	return std::move(multiply(m_1, m_2));
}

Matrix operator*(Matrix& m, double d)
{
	int numRows = m.getNumRows();
	int numColumns = m.getNumColumns();
	auto rose = makeIndexingSet(numRows);
	auto calls = makeIndexingSet(numColumns);

	auto md = makeMatrix_unique(numRows, numColumns);

	for (auto rho : rose)
		for (auto xi : calls)
			md[rho][xi] = m[rho][xi] * d;

	Matrix mtd(numRows, numColumns, std::move(md));

	return std::move(mtd);
}

Matrix operator*(double d, Matrix& m)
{
	int numRows = m.getNumRows();
	int numColumns = m.getNumColumns();
	auto rose = makeIndexingSet(numRows);
	auto calls = makeIndexingSet(numColumns);

	auto md = makeMatrix_unique(numRows, numColumns);

	for (auto rho : rose)
		for (auto xi : calls)
			md[rho][xi] = d * m[rho][xi];

	Matrix mtd(numRows, numColumns, std::move(md));

	return std::move(mtd);
}

Vector operator*(Matrix& A, Vector& x)
{
	int numRows = A.getNumRows();
	int numColumns = A.getNumColumns();
	assert(numColumns == x.getDimension());
	//if (numColumns != x.getDimension())
	//auto indices = makeIndexingSet(dim);
	auto rose = makeIndexingSet(numRows);
	auto calls = makeIndexingSet(numColumns);
	Vector b(numRows);

	for (auto rho : rose)
		for (auto xi : calls)
			b[rho] += A[rho][xi] * x[xi];

	return std::move(b);
}

bool operator==(Matrix& m_1, Matrix& m_2)
{
  return m_1.equals(m_2);
}

void Matrix::print()
{
	auto rowindices = makeIndexingSet(this->getNumRows());
	auto colindices = makeIndexingSet(this->getNumColumns());

  std::cout << "[";

  for (int i : rowindices)
  {
  	std::cout << "(";

  	for (int j : colindices)
  	{
    	if (j == (colindices.back() - 1))
      	std::cout << this->matrix[i][j];
    	else
      	std::cout << this->matrix[i][j] << ", ";
  	}

  	if (i == (rowindices.back() - 1))
  		std::cout << ")";
  	else
  		std::cout << ")," << std::endl;
  }

  std::cout << "]" << std::endl;
}

}
