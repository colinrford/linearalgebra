
#include "matrix.hpp"
#include <cstdint>
#include <ranges>
#include <iostream>

namespace linalg
{

auto make_indexing_set = [](int n) {
	return std::views::iota(0, n);
};

auto make_indexing_set_starting_at = [](int index, int size) {
	return std::views::iota(index, index + size);
};

auto make_indexing_set_reverse = [](int n) {
	return std::views::reverse(make_indexing_set(n));
};

auto make_matrix_unique = [](int numRows, int numColumns) {
	auto mtrx = std::make_unique<std::unique_ptr<double[]>[]>(numRows);
	auto rose = make_indexing_set(numRows);
	for (auto rho : rose)
		mtrx[rho] = std::make_unique<double[]>(numColumns);
	return mtrx;
};

/*
auto checkForExistingLUdcmp = [](Matrix& dis) {
	bool alreadyComputedLUdcmp = dis.lud.has_value();
	if (!alreadyComputedLUdcmp)
	{
		auto lu_decomp = this->croutLU();

		if (!lu_decomp.has_value())
			return alreadyComputedLUdcmp;

		LUdcmp decomp = std::move(lu_decomp.value());
		dis.lud = std::make_unique<LUdcmp>(std::move(decomp));
		alreadyComputedLUdcmp = true;
	}
	return alreadyComputedLUdcmp;
};
*/

/*constexpr auto zero_vector_of_dim = [](const std::size_t d) {
	vector v(d);
	for (const auto entry : v)
		entry = 0.;
	return v;
};*/

// The basic constructor just creates an n x m identity matrix
matrix::matrix(const std::size_t n, const std::size_t m)
{
	if (n == 0 || m == 0)
		throw matrix_exception::non_pos();

	n_rows = n;
	m_columns = m;
	const auto rose = make_indexing_set(n_rows);
	const auto calls = make_indexing_set(m_columns);

	matrix_ptr = make_matrix_unique(n, m);

	for (auto rho : rose)
		for (auto xi : calls)
			if (rho != xi)
				matrix_ptr[rho][xi] = 0.;
			else
				matrix_ptr[rho][xi] = 1.;
}

matrix::matrix(const std::size_t n, const std::size_t m,
								std::unique_ptr<std::unique_ptr<double[]>[]> mtrx)
								: matrix_ptr{std::move(mtrx)}
{
	if (n == 0 || m == 0)
		throw matrix_exception::non_pos();

	n_rows = n;
	m_columns = m;
}

matrix::matrix(const std::size_t n)
{
	if (n == 0)
		throw matrix_exception::non_pos();

	n_rows = m_columns = n;
	const auto square = make_indexing_set(n_rows);

	matrix_ptr = make_matrix_unique(n, n);

	for (auto rho : square)
		for (auto xi : square)
			if (rho != xi)
				matrix_ptr[rho][xi] = 0.;
			else
				matrix_ptr[rho][xi] = 1.;
}

matrix::matrix(const std::size_t n,
								std::unique_ptr<std::unique_ptr<double[]>[]> mtrx)
								: matrix_ptr{std::move(mtrx)}
{
	if (n == 0)
		throw matrix_exception::non_pos();

	n_rows = m_columns = n;
}

matrix::matrix(const std::vector<std::vector<double>> mtrx)
{
	if (mtrx.empty())
		throw matrix_exception::non_pos();
	n_rows = mtrx.size();
	m_columns = mtrx.front().size();
	const auto rose = make_indexing_set(n_rows);
	const auto calls = make_indexing_set(m_columns);

	matrix_ptr = make_matrix_unique(rose.size(), calls.size());

	for (auto rho : rose)
		for (auto xi : calls)
			matrix_ptr[rho][xi] = mtrx[rho][xi];
}

matrix::matrix(const std::vector<vector> rows)
{
	if (rows.empty())
		throw matrix_exception::non_pos();
	const std::size_t numColumns = rows.front().get_dimension();
	n_rows = rows.size();
	m_columns = numColumns;
	auto rose = make_indexing_set(n_rows);
	auto calls = make_indexing_set(m_columns);

	matrix_ptr = make_matrix_unique(rows.size(), calls.size());

	for (auto rho : rose)
		for (auto xi : calls)
			matrix_ptr[rho][xi] = rows[rho][xi];
}

matrix::matrix(const std::vector<double> diag)
{
	if (diag.empty())
		throw matrix_exception::non_pos();
	const std::size_t n = diag.size();
	n_rows = m_columns = n;
	const auto square = make_indexing_set(n);

	matrix_ptr = make_matrix_unique(n, n);

	for (auto rho : square)
		for (auto xi : square)
			if (rho != xi)
				matrix_ptr[rho][xi] = 0.;
			else
				matrix_ptr[rho][xi] = diag[rho];
}

matrix::matrix(const vector& diag)
{
	/*if (diag.empty())
		throw matrix_exception::non_pos();*/ // no check atm, assume nice input
	n_rows = m_columns = diag.size();
	const auto square = make_indexing_set(n_rows);

	matrix_ptr = make_matrix_unique(n_rows, n_rows);

	for (auto rho : square)
		for (auto xi : square)
			if (rho != xi)
				matrix_ptr[rho][xi] = 0.0;
			else
				matrix_ptr[rho][xi] = diag[rho];
}

matrix::matrix(matrix&& mtrx) noexcept
			: matrix_ptr{std::move_if_noexcept(mtrx.matrix_ptr)},
				n_rows{mtrx.n_rows},
				m_columns{mtrx.m_columns}
{
  mtrx.n_rows = mtrx.m_columns = 0; //unique_ptr.reset()?
}


matrix& matrix::operator=(matrix&& mtrx)
{
	this->matrix_ptr = std::move(mtrx.matrix_ptr);
	this->n_rows = mtrx.n_rows;
	this->m_columns = mtrx.m_columns;
	return *this;
}


auto matrix::operator[](const std::size_t index) -> decltype(matrix_ptr[index])
{
	if (index < this->get_num_rows())
		return this->matrix_ptr[index];
	else
		throw matrix_exception::out_of_bounds();
}

auto matrix::operator[](const std::size_t index) const
																								 -> decltype(matrix_ptr[index])
{
	if (index < this->get_num_rows())
		return this->matrix_ptr[index];
	else
		throw matrix_exception::out_of_bounds();
}

double& matrix::operator[](const std::size_t i, const std::size_t j)
{
	if (i < this->get_num_rows() && j < this->get_num_columns())
		return this->matrix_ptr[i][j];
	else
		throw matrix_exception::out_of_bounds();
}

double& matrix::operator[](const std::size_t i, const std::size_t j) const
{
	if (i < this->get_num_rows() && j < this->get_num_columns())
		return this->matrix_ptr[i][j];
	else
		throw matrix_exception::out_of_bounds();
}

matrix::~matrix() = default;

/*double** Matrix::getMatrix()
{
	//auto m = matrix.get();

	return this->matrix;
}*/
//suspicious
/*
void Matrix::setMatrix(std::unique_ptr<std::unique_ptr<double[]>[]> mtrx)
{
	matrix = std::move(mtrx);
}
*/

std::size_t matrix::get_num_rows() const { return this->n_rows; }

std::size_t matrix::get_num_columns() const { return this->m_columns; }

/*
void Matrix::setNumRows(std::size_t rows)
{
	this->n_rows = rows;
}

void Matrix::setNumColumns(std::size_t columns)
{
	this->m_columns = columns;
}
*/

double matrix::get_entry(const std::size_t i, const std::size_t j) const
{
	if ((i < this->get_num_rows()) && (j < this->get_num_columns()))
		return this->matrix_ptr[i][j];
	else
		throw matrix_exception::out_of_bounds();
}

/*
 *  iterator implementation
 */

struct matrix::row_iterator
{
public:
	using value_type = std::unique_ptr<double[]>;
	using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;
  using pointer = std::unique_ptr<double[]>*;
  using reference = std::unique_ptr<double[]>&;
  using iterator_category = std::random_access_iterator_tag;
	// does not actually compile with random_access std::algorithms as of 2/29/24
	// seems to qualify as, and work with, bidirectional std::algorithms

public:

	constexpr row_iterator() noexcept = default;

	constexpr row_iterator(pointer row_ptr, uint64_t ri) noexcept
	{
		iter_ptr = row_ptr;
		row_index = ri;
	}

	constexpr pointer operator->() noexcept { return iter_ptr; }

  constexpr reference operator*() noexcept { return *iter_ptr; }

  constexpr reference operator[](difference_type& index)
  { return *(iter_ptr + index); } // doesn't work

  constexpr row_iterator& operator++()
  {
		++row_index;
    ++iter_ptr;
    return *this;
  }

  constexpr row_iterator operator++(int)
  {
    auto iter = *this;
    ++(*this);
    return iter;
  }

  constexpr row_iterator& operator--()
  {
		--row_index;
    --iter_ptr;
    return *this;
  }

  constexpr row_iterator operator--(int)
  {
    auto iter = *this;
    --(*this);
    return iter;
  }

  constexpr row_iterator& operator+=(int off)
  {
    iter_ptr += off;
    return *this;
  }

  constexpr row_iterator operator+(difference_type& off)
  {
    auto iter = *this;
    return iter += off;// return iterator(iter + off);// based off gcc stl impl
  }

  constexpr row_iterator& operator-=(int off)
	{ return *this += -off; }

  constexpr row_iterator operator-(difference_type& off)
  {
    auto iter = *this;
    return iter -= off;
  }

  constexpr bool operator==(const row_iterator& other) const noexcept
  { return iter_ptr == other.iter_ptr; }

  constexpr bool operator!=(const row_iterator& other) const noexcept
  { return iter_ptr != other.iter_ptr; } // should it be !( == )?

  constexpr bool operator<(row_iterator other) const noexcept
  { return iter_ptr < other.iter_ptr; }

  constexpr bool operator>(row_iterator other) const noexcept
  { return iter_ptr > other.iter_ptr; }

  constexpr bool operator<=(row_iterator other) const noexcept
  { return iter_ptr <= other.iter_ptr; }

  constexpr bool operator>=(row_iterator other) const noexcept
  { return iter_ptr >= other.iter_ptr; }

  constexpr difference_type operator-(row_iterator other) const noexcept
  { return iter_ptr - other.iter_ptr; }

  constexpr friend auto operator<=>(row_iterator, row_iterator) = default;

	/*
  constexpr friend difference_type operator-(row_iterator first,
																						 row_iterator second)
  { return first - second; } // are

  constexpr friend row_iterator operator+(row_iterator it, int off)
  { return it += off; } // these

  constexpr friend row_iterator operator-(row_iterator it, int off)
  { return it -= off; } // okay? sure */
	friend struct row_major_iterator;
	friend struct column_major_iterator;

private:
	pointer iter_ptr;
	uint64_t row_index;
};

constexpr matrix::row_iterator matrix::row_begin()
{ return row_iterator(matrix_ptr.get(), 0); }

constexpr matrix::row_iterator matrix::row_end()
{ return row_iterator(matrix_ptr.get() + n_rows, n_rows); }

constexpr matrix::row_iterator row_begin(matrix& m)
{ return m.row_begin(); }

constexpr matrix::row_iterator row_end(matrix& m)
{ return m.row_end(); }


struct matrix::row_major_iterator
{
public:

  using value_type = double;
  using difference_type = std::ptrdiff_t;
  using pointer = double*;
  using reference = double&;
  using iterator_category = std::random_access_iterator_tag;

public:

  constexpr row_major_iterator() noexcept = default;

  constexpr row_major_iterator(pointer mtrx_ptr,
															 row_iterator rit,
														 	 matrix* which_matrix,
														 	 uint64_t ci) noexcept
							: current_row{rit},
								this_matrix{which_matrix}
	{
		iter_ptr = mtrx_ptr;
		col_index = ci;
	}

  constexpr pointer operator->() noexcept { return iter_ptr; }

  constexpr reference operator*() noexcept { return *iter_ptr; }

  constexpr reference operator[](difference_type index)
  {
		return *(iter_ptr + index);
	} // should this exist??

	//rough stuff
  constexpr row_major_iterator& operator++()
  {
		++col_index;
		if (col_index < this_matrix->m_columns)
    	++iter_ptr;
		else
		{
			++current_row;
			if (current_row.row_index < this_matrix->n_rows)
			{
				iter_ptr = current_row->get();
				col_index = 0;
			} else
				*this = this_matrix->row_major_end();
		}
    return *this;
  } //end rough stuff (haha just kidding)

  constexpr row_major_iterator operator++(int)
  {
    auto iter = *this;
    ++(*this);
    return iter;
  }

  constexpr row_major_iterator& operator--()
  {
		if (col_index != 0)
		{
			--col_index;
			--iter_ptr;
		}
		else
		{
			if (current_row.row_index != 0)
			{
				--current_row;
				col_index = this_matrix->m_columns - 1;
				iter_ptr = current_row->get() + col_index;
			} else
				*this = this_matrix->row_major_begin();
		}
		--current_row;
    return *this;
  } // alas

  constexpr row_major_iterator operator--(int)
  {
    auto iter = *this;
    --(*this);
    return iter;
  }

  constexpr row_major_iterator& operator+=(int off)
  {
    iter_ptr += off;
    return *this;
  }

  constexpr row_major_iterator operator+(difference_type& off)
  {
    auto iter = *this;
    return iter += off;// return iterator(iter + off);// based off gcc stl impl
  }

  constexpr row_major_iterator& operator-=(int off)
	{ return *this += -off; }

  constexpr row_major_iterator operator-(difference_type& off)
  {
    auto iter = *this;
    return iter -= off;
  }

  constexpr bool operator==(const row_major_iterator& other) const noexcept
  { return iter_ptr == other.iter_ptr; }

  constexpr bool operator!=(const row_major_iterator& other) const noexcept
  { return iter_ptr != other.iter_ptr; } // should it be !( == )?

  constexpr bool operator<(row_major_iterator other) const noexcept
  { return iter_ptr < other.iter_ptr; }

  constexpr bool operator>(row_major_iterator other) const noexcept
  { return iter_ptr > other.iter_ptr; }

  constexpr bool operator<=(row_major_iterator other) const noexcept
  { return iter_ptr <= other.iter_ptr; }

  constexpr bool operator>=(row_major_iterator other) const noexcept
  { return iter_ptr >= other.iter_ptr; }

  constexpr difference_type operator-(row_major_iterator other) const noexcept
  { return iter_ptr - other.iter_ptr; }

  constexpr friend auto operator<=>(row_major_iterator,
																		row_major_iterator) = default;

  constexpr friend difference_type operator-(row_major_iterator first,
																						 row_major_iterator second)
  { return *first - *second; } // are

  constexpr friend row_major_iterator operator+(row_major_iterator it, int off)
  { return it += off; } // these

  constexpr friend row_major_iterator operator-(row_major_iterator it, int off)
  { return it -= off; } // okay? sure

  //friend iterator operator+(int off, iterator);

private:
  pointer iter_ptr;
	matrix* this_matrix;
	row_iterator current_row;
	uint64_t col_index;
};

constexpr matrix::row_major_iterator matrix::row_major_begin()
{ return row_major_iterator(matrix_ptr.get()->get(), row_begin(), this, 0); }

constexpr matrix::row_major_iterator matrix::row_major_end()
{
	return row_major_iterator((matrix_ptr.get() + n_rows - 1)->get()
																											 + m_columns,
														 row_end(),
													 	 this,
													 	 this->m_columns);
}

matrix::row_major_iterator row_major_begin(matrix& m)
{ return m.row_major_begin(); }

matrix::row_major_iterator row_major_end(matrix& m)
{ return m.row_major_end(); }

// column_iterator is probably defunct
struct matrix::column_major_iterator
{
public:
	using value_type = double;
  using difference_type = std::ptrdiff_t;
  using pointer = double*;
  using reference = double&;
  using iterator_category = std::random_access_iterator_tag;

public:

	constexpr column_major_iterator() noexcept = default;

	constexpr column_major_iterator(pointer row_ptr,
																	matrix* which_matrix,
																	row_iterator row_tracker) noexcept
							: this_matrix{which_matrix},
								current_row{row_tracker}
	{
		iter_ptr = row_ptr;
		col_index = 0;
	}

	constexpr pointer operator->() noexcept { return iter_ptr; }

  constexpr reference operator*() noexcept { return *iter_ptr; }

  constexpr reference operator[](difference_type index)
  { return *(iter_ptr + index); }

  constexpr column_major_iterator& operator++()
  {
		++col_index;
		if (col_index < this_matrix->m_columns)
			++iter_ptr;
		else
			*this = this_matrix->column_major_end();
    return *this;
  }

  constexpr column_major_iterator operator++(int)
  {
    auto iter = *this;
    ++(*this);
    return iter;
  }

  constexpr column_major_iterator& operator--()
  {
		--col_index;
    --iter_ptr;
    return *this;
  }

  constexpr column_major_iterator operator--(int)
  {
    auto iter = *this;
    --(*this);
    return iter;
  }

  constexpr column_major_iterator& operator+=(int off)
  {
    iter_ptr += off;
    return *this;
  }

  constexpr column_major_iterator operator+(difference_type& off)
  {
    auto iter = *this;
    return iter += off;// return iterator(iter + off);// based off gcc stl impl
  }

  constexpr column_major_iterator& operator-=(int off)
	{ return *this += -off; }

  constexpr column_major_iterator operator-(difference_type& off)
  {
    auto iter = *this;
    return iter -= off;
  }

  constexpr bool operator==(const column_major_iterator& other) const noexcept
  { return iter_ptr == other.iter_ptr; }

  constexpr bool operator!=(const column_major_iterator& other) const noexcept
  { return iter_ptr != other.iter_ptr; } // should it be !( == )?

  constexpr bool operator<(column_major_iterator other) const noexcept
  { return iter_ptr < other.iter_ptr; }

  constexpr bool operator>(column_major_iterator other) const noexcept
  { return iter_ptr > other.iter_ptr; }

  constexpr bool operator<=(column_major_iterator other) const noexcept
  { return iter_ptr <= other.iter_ptr; }

  constexpr bool operator>=(column_major_iterator other) const noexcept
  { return iter_ptr >= other.iter_ptr; }

  constexpr difference_type
	operator-(column_major_iterator other) const noexcept
  { return iter_ptr - other.iter_ptr; }

  constexpr friend auto operator<=>(column_major_iterator,
																		column_major_iterator) = default;

	/*
  constexpr friend difference_type operator-(row_iterator first,
																						 row_iterator second)
  { return first - second; } // are

  constexpr friend row_iterator operator+(row_iterator it, int off)
  { return it += off; } // these

  constexpr friend row_iterator operator-(row_iterator it, int off)
  { return it -= off; } // okay? sure */
	//friend struct column_major_iterator;

private:
	pointer iter_ptr;
	matrix* this_matrix;
	row_iterator current_row;
	uint64_t col_index;
};

constexpr matrix::column_major_iterator matrix::column_major_begin()
{
	return column_major_iterator(matrix_ptr.get()->get(),
															 this,
															 this->row_begin());
}

constexpr matrix::column_major_iterator matrix::column_major_end()
{
	return column_major_iterator((matrix_ptr.get() + n_rows - 1)->get()
																								 + m_columns,
																this,
																this->row_end());
}

constexpr matrix::column_major_iterator column_major_begin(matrix& m)
{ return m.column_major_begin(); }

constexpr matrix::column_major_iterator column_major_end(matrix& m)
{ return m.column_major_end(); }

/* END ITERATOR IMPLEMENTATION */

bool can_add(const matrix& m_1, const matrix& m_2)
{
	const std::size_t rowsM_1 = m_1.get_num_rows();
	const std::size_t rowsM_2 = m_2.get_num_rows();

	const std::size_t columnsM_1 = m_1.get_num_columns();
	const std::size_t columnsM_2 = m_2.get_num_columns();

	if ((rowsM_1 != rowsM_2) || (columnsM_1 != columnsM_2))
		return false;
	else
		return true;
}

matrix matrix::add(const matrix& m_2) const
{
	if (!can_add(*this, m_2))
		throw matrix_exception::add_not_def();

	const std::size_t numRows = this->get_num_rows();
	const std::size_t numColumns = this->get_num_columns();
	const auto rose = make_indexing_set(numRows);
	const auto calls = make_indexing_set(numColumns);

	auto entries = make_matrix_unique(numRows, numColumns);

	for (auto rho : rose)
		for (auto xi : calls)
			entries[rho][xi] = this->matrix_ptr[rho][xi] + m_2.matrix_ptr[rho][xi];

	return matrix(numRows, numColumns, std::move(entries));
}

matrix add(const matrix& m_1, const matrix& m_2) { return m_1.add(m_2); }

bool can_subtract(const matrix& m_1, const matrix& m_2)
{ return can_add(m_1, m_2); }

matrix matrix::subtract(const matrix& m_2) const
{
	if (!can_subtract(*this, m_2))
		throw matrix_exception::sub_not_def();

	std::size_t numRows = this->get_num_rows();
	std::size_t numColumns = this->get_num_columns();
	auto rose = make_indexing_set(numRows);
	auto calls = make_indexing_set(numColumns);

	auto entries = make_matrix_unique(numRows, numColumns);

	for (auto rho : rose)
		for (auto xi : calls)
			entries[rho][xi] = this->matrix_ptr[rho][xi] - m_2.matrix_ptr[rho][xi];

	return matrix(numRows, numColumns, std::move(entries));
}

matrix subtract(const matrix& m_1, const matrix& m_2)
{ return m_1.subtract(m_2); }

bool can_multiply(const matrix& m_1, const matrix& m_2)
{
	const auto columnsM_1 = m_1.get_num_columns();
	const auto rowsM_2 = m_2.get_num_rows();

	if (columnsM_1 == rowsM_2)
		return true;
	else
		return false;
}

matrix matrix::multiply(const matrix& mp) const
{
	if (!can_multiply(*this, mp))
		throw matrix_exception::mult_not_def();

	const std::size_t numRowsNM = this->get_num_rows();
	const std::size_t numColumnsNM = this->get_num_columns();
	const std::size_t numRowsMP = mp.get_num_rows();
	assert(numColumnsNM == numRowsMP);
	const std::size_t numColumnsMP = mp.get_num_columns();

	const auto rowsNP = make_indexing_set(numRowsNM);
	const auto columnsNP = make_indexing_set(numColumnsMP);
	const auto columnsNM = make_indexing_set(numColumnsNM);

	auto entries = make_matrix_unique(rowsNP.size(), columnsNP.size());

	for (auto i : rowsNP)
		for (auto j : columnsNP)
			for (auto k : columnsNM)
				entries[i][j] += this->matrix_ptr[i][k] * mp[k, j];

	return matrix(rowsNP.size(), columnsNP.size(), std::move(entries));
}

matrix multiply(const matrix& nm, const matrix& mp)
{ return nm.multiply(mp); }

// incomplete
bool is_invertible(const matrix& nm)
{
	if (!nm.is_square())
		return false;
	/*if (compare(0, nm.det()))
		return false;*/

	return true;
}

// TODO: implement // incomplete
matrix matrix::inverse()
{
	if (!is_invertible(*this))
		throw matrix_exception::not_inv();

	const std::size_t numRows = this->get_num_rows();
	const std::size_t numColumns = this->get_num_columns();

	auto entries = make_matrix_unique(numRows, numColumns);

	return matrix(numRows, numColumns, std::move(entries));
}

// incomplete
void gaussje()
{

}

// incomplete
void gaussjInv()
{

}

// incomplete
matrix matrix::transpose()
{
	return matrix(this->get_num_columns(), this->get_num_rows());
}

// TODO: implement
double matrix::determinant()
{
	return 12345;
}

double matrix::det()
{
	return this->determinant();
}

// from Numerical Recipes, 3rd Ed – A is square
// std::optional w Matrix, optionally singular?
// Crout uses unit diagonals for the upper triangle U
bool matrix::croutLU()
{
	int imax = 0;
	int parity = 1;
	const std::size_t n = this->get_num_rows();
	double big;
	double temp;
	const auto square = make_indexing_set(n);
	vector perms(n);
	vector implicitScalingPerRow(n);
	auto lu = make_matrix_unique(n, n);

	for (auto i : square)
	{
		big = 0.;
		for (auto j : square)
		{
			lu[i][j] = this->matrix_ptr[i][j];
			if ((temp = abs(lu[i][j])) > big)
				big = temp;
		}
		if (compare(big, 0.))
		{
			this->lud = std::nullopt;
			return false; // matrix is singular
		}
		implicitScalingPerRow[i] = 1.0 / big;
	}
	for (auto k : square)
	{
		big = 0.;
		auto someIndices = make_indexing_set_starting_at(k, n - k);
		for (auto i : someIndices)
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
		perms[k] = imax;
		if (compare(lu[k][k], 0.)) // if 0, matrix is singular?
		{
			std::cout << "singular stuff goin on" << std::endl;
			lu[k][k] = std::numeric_limits<double>::epsilon();
		} // this if statement does not seem to behave quite right
		someIndices = make_indexing_set_starting_at(k + 1, n - (k + 1));
		for (auto i : someIndices)
		{
			temp = lu[i][k] /= lu[k][k];
			for (auto j : someIndices)
				lu[i][j] -= temp * lu[k][j];
		}
	}
	matrix lu_decomposition = matrix(n, std::move(lu));
	LUdcmp lu_decomp(std::move(lu_decomposition), std::move(perms), parity);
	this->lud = std::make_unique<LUdcmp>(std::move(lu_decomp));
	return true;
}

double matrix::croutLUDet()
{
	bool alreadyComputedLUdcmp = this->lud.has_value();
	if (!alreadyComputedLUdcmp)
		if (!this->croutLU())
			return 0.0;

	LUdcmp* decomp = this->lud.value().get();
	if (decomp->parity == 1 || decomp->parity == -1)
	{
		double dtrmnnt = decomp->parity;
		auto diag = make_indexing_set(this->n_rows);

		for (auto indx : diag)
			dtrmnnt *= decomp->lu_decomp[indx, indx];

		return dtrmnnt;
	}
	else
		return 0.0;
}

vector matrix::croutLUSolveSystem(const vector& b)
{
	bool alreadyComputedLUdcmp = this->lud.has_value();
	if (!alreadyComputedLUdcmp)
		if (!this->croutLU()) //std::optional?
			return vector(b.get_dimension());

	LUdcmp* decomp = this->lud.value().get();

	int ii = 0;
	int ip;
	int n = b.get_dimension();
	double sum = 0.0;
	vector x = make_zero_vector(n);

	matrix lu = 1 * decomp->lu_decomp;//!!!!!!!!!!!!!!!!!! copy ;0)
	auto vecIndices = make_indexing_set(n);

	for (auto rho : vecIndices)
		x[rho] = b[rho];

	for (auto i : vecIndices)
	{
		ip = decomp->permutations[i];
		sum = x[ip];
		x[ip] = x[i];
		if (ii != 0)
		{
			auto someIndices = make_indexing_set_starting_at(ii - 1, i - (ii - 1));
			for (auto j : someIndices)
				sum -= lu[i][j] * x[j];
		}
		else if (!compare(sum, 0.0))
    	ii = i + 1;
		x[i] = sum;
	}

	auto reversedIndices = make_indexing_set_reverse(n);
	for (auto i : reversedIndices)
	{
		sum = x[i];
		auto someIndices = make_indexing_set_starting_at(i + 1, n - (i + 1));
		for (auto j : someIndices)
			sum -= lu[i][j] * x[j];
		x[i] = sum / lu[i][i];
	}

	return x;
}

matrix matrix::croutLUSolveMatrixSystem(const matrix& B)
{
	const std::size_t numRowsA = this->get_num_rows();
	const int numColumnsA = this->get_num_columns();
	const int numRowsB = B.get_num_rows();
	const int numColumnsB = B.get_num_columns();

	if (numRowsA != numRowsB || numColumnsA != numColumnsB)
	{ // std::optional instead?
		//throw("LUdcmp::solve bad sizes");
		return matrix(numRowsA, numColumnsA);
	}

	auto rose = make_indexing_set(numRowsA);
	auto calls = make_indexing_set(numColumnsA);

	matrix X(numRowsA, numColumnsA);
	vector b = make_zero_vector(numRowsA);
	vector x = make_zero_vector(numRowsA);
	for (auto j : calls)
	{
		for (auto rho : rose)
			b[rho] = B[rho][j];
		x = std::move(croutLUSolveSystem(b));
		for (auto rho : rose)
			X[rho, j] = x[rho];
	}

	return X;
}

matrix matrix::croutLUInv()
{
	bool alreadyComputedLUdcmp = this->lud.has_value();
	if (!alreadyComputedLUdcmp)
		if (!this->croutLU()) //std::optional?
			return matrix(this->get_num_rows(), this->get_num_columns());

	LUdcmp* decomp = this->lud.value().get();

	std::size_t n = decomp->lu_decomp.get_num_rows();
	auto square = make_indexing_set(n);
	auto ainv_entries = make_matrix_unique(n, n);

  for (auto i : square)
		for (auto j : square)
			if (i != j)
				ainv_entries[i][j] = 0.0;
			else
				ainv_entries[i][j] = 1.0;

	matrix ainv(n, std::move(ainv_entries));
  matrix A_inv = this->croutLUSolveMatrixSystem(ainv);
	return A_inv;
}

// adapted from:	sci.utah.edu/~wallstedt
// main diagonal of L is composed with 1s
std::optional<LUdcmp> matrix::doolittleLU()
{
	const std::size_t n = this->get_num_rows();
	const auto square = make_indexing_set(n);
	auto lu = make_matrix_unique(n, n);

	for (auto k : square)
	{
		const auto subsq = make_indexing_set(k + 1);
		const auto subsqc = make_indexing_set_starting_at(k, n - k);
		const auto subsqcpo = make_indexing_set_starting_at(k + 1, n - (k + 1));
    for (auto j : subsqc)//(int j = k; j < n; ++j)
		{
      double sum = 0.0;
      for (auto p : subsq)
				sum += lu[k][p] * lu[p][j];
      lu[k][j] = this->matrix_ptr[k][j] - sum;
    }
    for (auto i : subsqcpo)//(int i = k + 1; i < n; ++i)
		{
    	double sum = 0.0;
    	for (auto p : subsq)
				sum += lu[i][p] * lu[p][k];
    	lu[i][k] = (this->matrix_ptr[i][k] - sum) / lu[k][k];
    }
	}
	matrix lu_decomposition(n, std::move(lu));
	vector perms = make_zero_vector(n);
	int parity = 0;

	return LUdcmp(std::move(lu_decomposition), std::move(perms), parity);
}

//incomplete possibly may need permutation info too?
double matrix::doolittleLUDet()
{
	auto lu_decomp = this->doolittleLU();

	if (!lu_decomp.has_value()) //std::optional?
		return 0.0;

	LUdcmp decomp = std::move(lu_decomp.value());

	double dtrmnnt = 1.0;
	auto diag = make_indexing_set(this->n_rows);

	for (auto indx : diag)
		dtrmnnt *= decomp.lu_decomp[indx][indx];

	return dtrmnnt;
}

bool matrix::is_square() const
{
	if (this->get_num_rows() == this->get_num_columns())
		return true;
	else
		return false;
}

bool matrix::equals(const matrix& m_2) const
{
	const std::size_t numRows_M_1 = this->get_num_rows();
	const std::size_t numColumns_M_1 = this->get_num_columns();
	const std::size_t numRows_M_2 = m_2.get_num_rows();
	const std::size_t numColumns_M_2 = m_2.get_num_columns();

  if ((numRows_M_1 != numRows_M_2) || (numColumns_M_1 != numColumns_M_2))
    return false;

	auto rose = make_indexing_set(numRows_M_1);
	auto calls = make_indexing_set(numColumns_M_1);

  for (auto rho : rose)
  	for (auto xi : calls)
    	if (!compare(this->matrix_ptr[rho][xi], m_2[rho, xi]))
     		return false;

  return true;
}

matrix operator+(const matrix& m_1, const matrix& m_2)
{ return add(m_1, m_2); }

matrix operator-(const matrix& m_1, const matrix& m_2)
{ return subtract(m_1, m_2); }

matrix operator*(const matrix& m_1, const matrix& m_2)
{ return multiply(m_1, m_2); }

matrix operator*(const matrix& m, const double d)
{
	const std::size_t numRows = m.get_num_rows();
	const std::size_t numColumns = m.get_num_columns();
	const auto rose = make_indexing_set(numRows);
	const auto calls = make_indexing_set(numColumns);

	auto md = matrix(numRows, numColumns);

	for (auto rho : rose)
		for (auto xi : calls)
			md[rho, xi] = m[rho, xi] * d;

	return md;
}

matrix operator*(const double d, const matrix& m)
{
	const std::size_t numRows = m.get_num_rows();
	const std::size_t numColumns = m.get_num_columns();
	const auto rose = make_indexing_set(numRows);
	const auto calls = make_indexing_set(numColumns);

	auto md = matrix(numRows, numColumns);
	for (auto rho : rose)
		for (auto xi : calls)
			md[rho, xi] = d * m[rho, xi];

	return md;
}

vector operator*(const vector& x, const matrix& A)
{
	const std::size_t numRows = A.get_num_rows();
	const std::size_t numColumns = A.get_num_columns();
	//if (numRows != x.get_dimension())
	//auto indices = make_indexing_set(dim);
	const auto rose = make_indexing_set(numRows);
	const auto calls = make_indexing_set(numColumns);
	vector b = make_zero_vector(numRows);

	for (auto rho : rose)
		for (auto xi : calls)
			b[rho] += x[rho] * A[rho, xi];

	return b;
}

vector operator*(const matrix& A, const vector& x)
{
	const std::size_t numRows = A.get_num_rows();
	const std::size_t numColumns = A.get_num_columns();
	assert(numColumns == x.get_dimension());
	//if (numColumns != x.get_dimension())
	//auto indices = make_indexing_set(dim);
	const auto rose = make_indexing_set(numRows);
	const auto calls = make_indexing_set(numColumns);
	vector b = make_zero_vector(numRows);

	for (auto rho : rose)
		for (auto xi : calls)
			b[rho] += A[rho, xi] * x[xi];

	return b;
}

bool operator==(const matrix& m_1, const matrix& m_2)
{ return m_1.equals(m_2); }

void matrix::print() const
{
	auto rowindices = make_indexing_set(this->get_num_rows());
	auto colindices = make_indexing_set(this->get_num_columns());

  std::cout << "[";

  for (int i : rowindices)
  {
  	std::cout << "(";

  	for (int j : colindices)
    	if (j != colindices.back())
      	std::cout << this->matrix_ptr[i][j] << ", ";
    	else
      	std::cout << this->matrix_ptr[i][j];

  	if (i != rowindices.back())
  		std::cout << ")," << std::endl;
  	else
  		std::cout << ")";
  }

  std::cout << "]" << std::endl;
}

//not written atm
void matrix::writeTeXto(std::string feyell)
{
	auto rowindices = make_indexing_set(this->get_num_rows());
	auto colindices = make_indexing_set(this->get_num_columns());

  std::cout << "[";

  for (int i : rowindices)
  {
  	std::cout << "(";

  	for (int j : colindices)
    	if (j != colindices.back())
      	std::cout << this->matrix_ptr[i][j] << ", ";
    	else
      	std::cout << this->matrix_ptr[i][j];

  	if (i != rowindices.back())
  		std::cout << ")," << std::endl;
  	else
  		std::cout << ")";
  }

  std::cout << "]" << std::endl;
}

}
