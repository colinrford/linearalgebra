
#include "matrix.h"

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

	std::vector<std::vector<double> > matrix(nRows, std::vector<double>(mColumns));

}

// Create matrix from given matrix data
Matrix::Matrix(int n, int m, std::vector<std::vector<double>> mtrx) : matrix{std::move(mtrx)}
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


auto Matrix::operator[](int index) -> decltype(this->getMatrix()[index])
{
	int numColumns = this->getNumColumns();
	if (index >= numColumns)
		throw MatrixException();

	return this->getMatrix()[index];
}

// I'd like this to be a (non-owning) raw pointer, if possible.
// Anything acceptable that is lightweight!!
std::vector<std::vector<double> > Matrix::getMatrix()
{
	std::vector<std::vector<double> > m = matrix;

	return m;
}

// 
void Matrix::setMatrix(std::vector<std::vector<double> > mtrx)
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

	std::vector<double> row = this->getMatrix()[i];
	return row[j];
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

// 
Matrix Matrix::add(Matrix& m_2)
{
	if (!canAdd(*this, m_2))
		throw MatrixException();

	int numRows = this->getNumRows();
	int numColumns = this->getNumColumns();

	std::vector<std::vector<double> > entries(numRows, std::vector<double>(numColumns));

	for (int i = 0; i < numRows; i++)
		for (int j = 0; j < numColumns; j++)
			entries[i][j] = this->getEntry(i, j) + m_2.getEntry(i, j);

	Matrix m1pm2(numRows, numColumns, std::move(entries));

	return m1pm2;
}

Matrix add(Matrix& m_1, Matrix& m_2)
{
	return m_1.add(m_2);
}

bool canSubtract(Matrix& m_1, Matrix& m_2)
{
	return canAdd(m_1, m_2);
}

// 
Matrix Matrix::subtract(Matrix& m_2)
{
	if (!canSubtract(*this, m_2))
		throw MatrixException();

	int numRows = this->getNumRows();
	int numColumns = this->getNumColumns();

	std::vector<std::vector<double> > entries(numRows, std::vector<double>(numColumns));

	for (int i = 0; i < numRows; i++)
		for (int j = 0; j < numColumns; j++)
			entries[i][j] = this->getEntry(i, j) - m_2.getEntry(i, j);

	Matrix m1pm2(numRows, numColumns, std::move(entries));

	return m1pm2;
}

Matrix subtract(Matrix& m_1, Matrix& m_2)
{
	return m_1.subtract(m_2);
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

// Matrix multiplication
Matrix Matrix::multiply(Matrix& mp)
{
	if (!canMultiply(*this, mp))
		throw MatrixException();

	int numRowsNM = this->getNumRows();
	int numColumnsNM = this->getNumColumns();
	int numRowsMP = mp.getNumRows();
	assert(numColumnsNM == numRowsMP);
	int numColumnsMP = mp.getNumColumns();

	std::vector<std::vector<double> > nm = this->getMatrix();
	std::vector<std::vector<double> > _mp = mp.getMatrix();
	std::vector<std::vector<double> > entries(numRowsNM, std::vector<double>(numColumnsMP));

	/* placeholder lol 
	for (this : that) 
	{
		for (that : this)
		{
	
		}
	}
	*/

	for (int i = 0; i < numRowsNM; i++)
	{
		for (int j = 0; j < numColumnsMP; j++)
		{
			for (int k = 0; k < numColumnsNM; k++)
			{
				entries[i][j] += nm[i][k] * _mp[k][j];
			}
		}
	}

	Matrix np(numRowsNM, numColumnsMP, std::move(entries));

	return np; 
}

Matrix multiply(Matrix& nm, Matrix& mp)
{
	return nm.multiply(mp);
}

bool isInvertible(Matrix& nm)
{
	if (!nm.isSquare())
	{
		throw MatrixException();
		return false;
	}
	if (compare(0, nm.det()))
	{
		throw MatrixException();
		return false;
	}

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

	std::vector<std::vector<double> > entries(numRows, std::vector<double>(numColumns));
	
	Matrix m(numRows, numColumns);

	return m;
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
 *        The permutation matrix is not stored as a matrix, but in an integer vector P of size N+1 
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

// Check if the given matrix is square
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

  for (int i = 0; i < numRowsM_1; i++)
  {
  	for (int j = 0; j < numColumnsM_1; j++)
  	{
    	if (this->matrix[i][j] != m_2.matrix[i][j])
    	{
     		return false;
    	}
    }
  }
  
  return true;
}

Matrix operator+(Matrix& m_1, Matrix& m_2)
{
  return add(m_1, m_2);
}

Matrix operator-(Matrix& m_1, Matrix& m_2)
{
  return subtract(m_1, m_2);
}

Matrix operator*(Matrix& m_1, Matrix& m_2)
{
	return multiply(m_1, m_2);
}

bool operator==(Matrix& m_1, Matrix& m_2)
{
  return m_1.equals(m_2);
}

void Matrix::print()
{
  std::vector<std::vector<double> > mtrx = this->getMatrix();

  int numRows = this->getNumRows();
  int numColumns = this->getNumColumns();

  std::cout << "[";

  // i has expected type of std::vector<T>. Think of it as one of the rows in mtrx
  for (int i = 0; i < numRows; i++)
  {
  	std::cout << "(";

  	// j has expected type of T (e.g. double). An entry of mtrx
  	for (int j = 0; j < numColumns; j++)
  	{
    	if (j == (numColumns - 1))
      	std::cout << mtrx[i][j];
    	else
      	std::cout << mtrx[i][j] << ", ";
  	}

  	if (i == (numRows - 1))
  		std::cout << ")";
  	else
  		std::cout << "),\n";
  }

  std::cout << "]\n";
}

bool compare(double a, double b)
{
  double epsilon = 1E-40;

  if (abs(b - a) < epsilon)
    return true;
  else
    return false;
}