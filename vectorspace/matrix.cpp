
#include "matrix.h"
#include <cassert>
#include <list>
#include <numeric>

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

Matrix::Matrix(std::vector<std::vector<double> > mtrx)
{
	auto rose = makeIndexingSet(mtrx.size());
	auto calls = makeIndexingSet(mtrx.front().size());

	matrix = makeMatrix_unique(rose.size(), calls.size());

	for (auto rho : rose)
		for (auto xi : calls)
			matrix[rho][xi] = mtrx[rho][xi];
}
/*
Matrix::Matrix(std::vector<Vector> rows)
{
	auto firstRow = &rows.front();
	int numColumns = firstRow->size();
	auto rose = makeIndexingSet(rows.size());
	auto calls = makeIndexingSet(numColumns);

	matrix = std::make_unique<std::unique_ptr<double[]>[]>(rows.size());
	for (auto rho : rose)
		matrix[rho] = std::make_unique<double[]>(calls.size());

	for (auto rho : rose)
		for (auto xi : calls)
			matrix[rho][xi] = rows[rho][xi];
}
*/
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
			entries[rho][xi] = this->getEntry(rho, xi) + m_2.getEntry(rho, xi);

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
			entries[rho][xi] = this->getEntry(rho, xi) - m_2.getEntry(rho, xi);

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

// TODO: implement
double Matrix::determinant()
{
	return 12345;
}

double Matrix::det()
{
	return this->determinant();
}


// LUP functions below based on C  code from Wikipedia page /wiki/LU_decomposition#Algorithms

/* INPUT: A - array of pointers to rows of a square matrix having dimension N
 *        Tol - small tolerance number to detect failure when the matrix is near degenerate
 * OUTPUT: Matrix A is changed, it contains a copy of both matrices L-E and U as A=(L-E)+U such that P*A=L*U.
 *        The permutation matrix is not stored as a matrix, but in an integer unique_ptr P of size N+1
 *        containing column indexes where the permutation matrix has "1". The last element P[N]=S+N,
 *        where S is the number of row exchanges needed for determinant computation, det(P)=(-1)^S
 */

/*
int LUPDecompose(double **A, int N, double Tol, int *P)
{

    int i, j, k, imax;
    double maxA, *ptr, absA;

    for (i = 0; i <= N; i++)
    {
        P[i] = i; // Unit permutation matrix, P[N] initialized with N
    }

    for (i = 0; i < N; i++)
    {
        maxA = 0.0;
        imax = i;

        for (k = i; k < N; k++)
        {
           if ((absA = fabs(A[k][i])) > maxA)
           {
           		maxA = absA;
           		imax = k;
           }
         }

        if (maxA < Tol)
        {
					return 0; // failure, matrix is degenerate
				}

        if (imax != i)
        {
            // pivoting P
            j = P[i];
            P[i] = P[imax];
            P[imax] = j;

            // pivoting rows of A
            ptr = A[i];
            A[i] = A[imax];
            A[imax] = ptr;

            // counting pivots starting from N (for determinant)
            P[N]++;
        }

        for (j = i + 1; j < N; j++)
        {
        	A[j][i] /= A[i][i];

          for (k = i + 1; k < N; k++)
          {
          	A[j][k] -= A[j][i] * A[i][k];
          }
        }
    }

    return 1;  // decomposition done
}
*/

/* INPUT: A,P filled in LUPDecompose; N - dimension.
 * OUTPUT: Function returns the determinant of the initial matrix
 */
/*
double LUPDeterminant(Matrix& A, int* P, int N)
{

	double det = A[0][0];

  for (int i = 1; i < N; i++)
  	det *= A[i][i];

  if ((P[N] - N) % 2 == 0)
    return det;
  else
    return -det;
}
*/

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

bool compare(double a, double b)
{
  double epsilon = 2 * std::numeric_limits<double>::epsilon();

  if (abs(b - a) < epsilon)
    return true;
  else
    return false;
}
