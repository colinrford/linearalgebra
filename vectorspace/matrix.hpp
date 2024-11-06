
#include "vector.hpp"
#include <cassert>
#include <list>
#include <numeric>
#include <string>
#include <utility>

namespace linalg
{

struct LUdcmp;

// also stale
struct matrix_exception : public std::exception {

	struct non_pos : public std::exception {
		const char* what () const throw () {
  		return "The number of rows and columns of the matrix must be at least 1";
  	}
	};

	struct add_not_def : public std::exception {
	  const char* what() const throw () {
	  	return "Addition is undefined for these matrices. The numbers of rows \
			 				and columns must be equal, i.e., n x m + n x m";
	  }
	};

	struct sub_not_def : public std::exception {
	  const char* what() const throw () {
	  	return "Subtraction is undefined for these matrices. The numbers of \
							rows and columns must be equal, i.e., n x m - n x m";
	  }
	};

  struct mult_not_def : public std::exception {
		const char* what() const throw () {
  		return "Multiplication is undefined for these matrices. Matrices must be \
						of the form n x m * m x p";
  	}
	};

	struct not_inv : public std::exception {
		const char* what() const throw () {
	  	return "This matrix is not invertible. Check: Is it square? \
							Is its determinant nonzero?";
	  }
	};

	struct out_of_bounds : public std::exception {
    const char* what () const throw () {
      return "Index exceeds dimension of matrix.\n";
    }
  };
};


class matrix {

	private:

		std::size_t n_rows;
		std::size_t m_columns;
		std::optional<std::string_view> label;
		std::unique_ptr<std::unique_ptr<double[]>[]> matrix_ptr;
		std::optional<std::unique_ptr<LUdcmp>> lud;

		matrix(const std::size_t n, const std::size_t m,
						std::unique_ptr<std::unique_ptr<double[]>[]> mtrx);
		matrix(const std::size_t n,
						std::unique_ptr<std::unique_ptr<double[]>[]> mtrx);



	public:

		matrix(const std::size_t n, const std::size_t m);

		matrix(const std::size_t n);

		matrix(const vector& diag);

		matrix(const std::vector<double> diag);

		matrix(const std::vector<std::vector<double>> mtrx);

		matrix(std::vector<linalg::vector> rows);

		matrix(matrix&& mtrx) noexcept;

		matrix& operator=(matrix&& mtrx);

		auto operator[](const std::size_t index) -> decltype(matrix_ptr[index]);

		auto operator[](const std::size_t index) const -> decltype(matrix_ptr[index]);

		double& operator[](const std::size_t i, const std::size_t j);

		double& operator[](const std::size_t i, const std::size_t j) const;

		std::size_t get_num_rows() const;

		std::size_t get_num_columns() const;

		double get_entry(const std::size_t i, const std::size_t j) const;

		struct row_iterator;
		constexpr row_iterator row_begin();
		constexpr row_iterator row_end();

		struct row_major_iterator;
		constexpr row_major_iterator row_major_begin();
		constexpr row_major_iterator row_major_end();

		struct column_iterator;
		constexpr column_iterator column_begin();
		constexpr column_iterator column_end();
		constexpr column_iterator col_begin();
		constexpr column_iterator col_end();

		struct column_major_iterator;
		constexpr column_major_iterator column_major_begin();
		constexpr column_major_iterator column_major_end();
		constexpr column_major_iterator col_major_begin();
		constexpr column_major_iterator col_major_end();

		//double** getmatrix();

		matrix add(const matrix& m_2) const;

		matrix subtract(const matrix& m_2) const;

		matrix multiply(const matrix& mp) const;

		matrix inverse();

		matrix transpose();

 		bool croutLU();

		vector croutLUSolveSystem(const vector&);

		matrix croutLUSolveMatrixSystem(const matrix&);

		double croutLUDet();

		matrix croutLUInv();

		std::optional<linalg::LUdcmp> doolittleLU();

		double doolittleLUDet();

		double determinant();

		double det();

		bool is_square() const;

		bool equals(const matrix&) const;

		void print() const;

		void writeTeXto(std::string);

		~matrix();// = default; // does this not work?

		struct iterator;

    iterator begin();

    iterator end();

};

struct LUdcmp
{
	matrix lu_decomp;
	vector permutations;
	int parity;

	LUdcmp(matrix lud, vector perms, int par)
					: lu_decomp(std::move(lud)),
						permutations(std::move(perms)),
						parity(par) { };

	double operator[](const std::size_t i, const std::size_t j)
	{ return lu_decomp[i, j]; };

	double operator[](const std::size_t i, const std::size_t j) const
	{ return lu_decomp[i, j]; };
};

struct matrix_info
{
	double determinant;
	matrix inverse;
	matrix transpose;
};

matrix operator+(const matrix&, const matrix&);
matrix operator-(const matrix&, const matrix&);
matrix operator*(const matrix&, const matrix&);
matrix operator*(const double, const matrix&);
matrix operator*(const matrix&, const double);
vector operator*(const matrix&, const vector&);
matrix operator/(const matrix&, const double);
bool operator==(const matrix&, const matrix&);
//void print(const matrix&);

}
