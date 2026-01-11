import std;
import lam.linearalgebra;
import lam.concepts;

using namespace lam::linalg;
using namespace lam::linalg::concepts::experimental;

// Mock "Foreign" Matrix (e.g., behaving like Eigen or another library)
// Uses rows()/cols() instead of get_rows()/get_cols()
// Uses operator() instead of operator[]
struct foreign_matrix
{
  using scalar_type = double;
  using value_type = double;
  using size_type = std::size_t;

  std::vector<double> data;
  size_type r, c;

  foreign_matrix(size_type rows, size_type cols) : r(rows), c(cols), data(rows * cols) {}

  // Required for additive_group_element_c_weak (T(0))
  // Represents a zero element (here, a 0x0 matrix for simplicity, or could be context dependent)
  explicit foreign_matrix(int) : r(0), c(0), data() {}

  // "Foreign" API
  size_type rows() const { return r; }
  size_type cols() const { return c; }
  double& operator()(size_type i, size_type j) { return data[i * c + j]; }
  const double& operator()(size_type i, size_type j) const { return data[i * c + j]; }

  // Add C++23 style subscript operator to fully satisfy matrix_c_weak
  double& operator[](size_type i, size_type j) { return data[i * c + j]; }
  const double& operator[](size_type i, size_type j) const { return data[i * c + j]; }

  // Add row/col accessors to fully satisfy matrix_c_weak
  std::span<double> row(size_type i) { return std::span<double>(data.data() + i * c, c); }
  auto col(size_type j) { return std::views::iota(0, 0); } // Dummy implementation

  // Operators to satisfy vectorspace_element_c_weak (Group + Module ops)
  foreign_matrix& operator+=(const foreign_matrix&) { return *this; }
  friend foreign_matrix operator+(foreign_matrix lhs, const foreign_matrix& rhs) { return lhs; }

  foreign_matrix& operator-=(const foreign_matrix&) { return *this; }
  friend foreign_matrix operator-(foreign_matrix lhs, const foreign_matrix& rhs) { return lhs; }

  foreign_matrix operator-() const { return *this; }

  friend foreign_matrix operator*(double, foreign_matrix m) { return m; }
  friend foreign_matrix operator*(foreign_matrix m, double) { return m; }
};

// 1. Assert it NOW PASSES the matrix concept because the API (rows/cols) matches!
// We have relaxed the exclusivity by adopting standard names.
static_assert(matrix_c_weak<foreign_matrix, double>,
              "foreign_matrix SHOULD satisfy matrix_c_weak now that we use rows()/cols()");

// 2. Assert it passes the generic vectorspace concept
static_assert(vectorspace_element_c_weak<foreign_matrix, double>,
              "foreign_matrix SHOULD satisfy vectorspace_element_c_weak");

int main()
{
  std::println("Exclusivity check updated:");
  std::println(" - foreign_matrix IS now a matrix_c_weak (standard API adoption)");
  std::println(" - foreign_matrix IS a vectorspace_element_c_weak");
  return 0;
}
