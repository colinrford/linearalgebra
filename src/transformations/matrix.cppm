export module lam.linearalgebra:matrix;
import std;
import lam.concepts;
import :vectorspace;
import :transformations.concepts;

namespace lam::linalg
{

export struct matrix_exception : public std::exception
{
  struct non_pos : public std::exception
  {
    const char* what() const noexcept override { return "Matrix dimensions must be at least 1x1."; }
  };
  struct dim_mismatch : public std::exception
  {
    const char* what() const noexcept override { return "Matrix dimensions mismatch for the operation."; }
  };
  struct out_of_bounds : public std::exception
  {
    const char* what() const noexcept override { return "Index out of bounds."; }
  };
  struct not_invertible : public std::exception
  {
    const char* what() const noexcept override { return "Matrix is not invertible."; }
  };
};

export template<typename T, typename Alloc = std::allocator<T>>
  requires lam::concepts::experimental::ring_element_c_weak<T>
class matrix
{
public:
  using value_type = T;
  using scalar_type = T; // Added for concept compatibility
  using allocator_type = Alloc;
  using size_type = std::size_t;

private:
  size_type rows;
  size_type cols;
  std::unique_ptr<T[], std::function<void(T*)>> data;
  std::optional<std::string> label;

  struct Deleter
  {
    Alloc alloc;
    size_type n;
    void operator()(T* p)
    {
      if (p)
        std::allocator_traits<Alloc>::deallocate(alloc, p, n);
    }
  };

public:
  constexpr matrix() : rows{0}, cols{0}, data{nullptr, [](T*) {}} {}

  matrix(size_type r, size_type c, const Alloc& alloc = Alloc())
    : rows{r}, cols{c}, data{nullptr, Deleter{alloc, r * c}}
  {
    if (r == 0 || c == 0)
      throw matrix_exception::non_pos();
    size_type n = r * c;
    auto* ptr = std::allocator_traits<Alloc>::allocate(const_cast<Alloc&>(alloc), n);
    data = std::unique_ptr<T[], std::function<void(T*)>>(ptr, Deleter{alloc, n});
    std::uninitialized_fill_n(ptr, n, T{0});
  }

  explicit matrix(size_type n, const Alloc& alloc = Alloc()) : matrix(n, n, alloc)
  {
    for (size_type i = 0; i < n; ++i)
      (*this)[i, i] = T{1};
  }

  matrix(const std::vector<std::vector<T>>& mtrx, const Alloc& alloc = Alloc())
    : matrix(mtrx.size(), mtrx.empty() ? 0 : mtrx[0].size(), alloc)
  {
    for (size_type i = 0; i < rows; ++i)
    {
      for (size_type j = 0; j < cols; ++j)
      {
        (*this)[i, j] = mtrx[i][j];
      }
    }
  }

  matrix(matrix&&) noexcept = default;
  matrix& operator=(matrix&&) noexcept = default;

  matrix(const matrix& other)
    : rows{other.rows}, cols{other.cols}, data{nullptr, Deleter{Alloc(), other.rows * other.cols}}
  {
    size_type n = rows * cols;
    auto* ptr = std::allocator_traits<Alloc>::allocate(Deleter{Alloc(), n}.alloc, n);
    data = std::unique_ptr<T[], std::function<void(T*)>>(ptr, Deleter{Alloc(), n});
    std::uninitialized_copy(other.data.get(), other.data.get() + n, ptr);
  }

  matrix& operator=(const matrix& other)
  {
    if (this != &other)
    {
      matrix tmp(other);
      *this = std::move(tmp);
    }
    return *this;
  }

  // C++23 Multidimensional []
  T& operator[](size_type i, size_type j)
  {
    if (i >= rows || j >= cols)
      throw matrix_exception::out_of_bounds();
    return data[i * cols + j];
  }

  const T& operator[](size_type i, size_type j) const
  {
    if (i >= rows || j >= cols)
      throw matrix_exception::out_of_bounds();
    return data[i * cols + j];
  }

  constexpr size_type get_rows() const noexcept { return rows; }
  constexpr size_type get_cols() const noexcept { return cols; }

  matrix transpose() const
  {
    matrix res(cols, rows);
    for (size_type i = 0; i < rows; ++i)
    {
      for (size_type j = 0; j < cols; ++j)
      {
        res[j, i] = (*this)[i, j];
      }
    }
    return res;
  }

  // ========== Iterators and Accessors ==========

  // Element iteration (row-major order)
  T* begin() noexcept { return data.get(); }
  T* end() noexcept { return data.get() + rows * cols; }
  const T* begin() const noexcept { return data.get(); }
  const T* end() const noexcept { return data.get() + rows * cols; }
  const T* cbegin() const noexcept { return data.get(); }
  const T* cend() const noexcept { return data.get() + rows * cols; }

  // Single row access — returns contiguous span
  std::span<T> row(size_type i)
  {
    if (i >= rows)
      throw matrix_exception::out_of_bounds();
    return std::span<T>{data.get() + i * cols, cols};
  }

  std::span<const T> row(size_type i) const
  {
    if (i >= rows)
      throw matrix_exception::out_of_bounds();
    return std::span<const T>{data.get() + i * cols, cols};
  }

  // All rows — returns range of spans
  auto rows_range()
  {
    return std::views::iota(size_type{0}, rows)
         | std::views::transform([this](size_type i) {
             return std::span<T>{data.get() + i * cols, cols};
           });
  }

  auto rows_range() const
  {
    return std::views::iota(size_type{0}, rows)
         | std::views::transform([this](size_type i) {
             return std::span<const T>{data.get() + i * cols, cols};
           });
  }

  // Single column access — strided view (non-contiguous)
  auto col(size_type j)
  {
    if (j >= cols)
      throw matrix_exception::out_of_bounds();
    return std::views::iota(size_type{0}, rows)
         | std::views::transform([this, j](size_type i) -> T& {
             return data[i * cols + j];
           });
  }

  auto col(size_type j) const
  {
    if (j >= cols)
      throw matrix_exception::out_of_bounds();
    return std::views::iota(size_type{0}, rows)
         | std::views::transform([this, j](size_type i) -> const T& {
             return data[i * cols + j];
           });
  }

  // All columns — returns range of strided ranges
  auto cols_range()
  {
    return std::views::iota(size_type{0}, cols)
         | std::views::transform([this](size_type j) {
             return std::views::iota(size_type{0}, rows)
                  | std::views::transform([this, j](size_type i) -> T& {
                      return data[i * cols + j];
                    });
           });
  }

  auto cols_range() const
  {
    return std::views::iota(size_type{0}, cols)
         | std::views::transform([this](size_type j) {
             return std::views::iota(size_type{0}, rows)
                  | std::views::transform([this, j](size_type i) -> const T& {
                      return data[i * cols + j];
                    });
           });
  }

  matrix& operator+=(const matrix& other)
  {
    if (rows != other.rows || cols != other.cols)
      throw matrix_exception::dim_mismatch();
    size_type n = rows * cols;
    for (size_type i = 0; i < n; ++i)
      data[i] += other.data[i];
    return *this;
  }

  matrix& operator-=(const matrix& other)
  {
    if (rows != other.rows || cols != other.cols)
      throw matrix_exception::dim_mismatch();
    size_type n = rows * cols;
    for (size_type i = 0; i < n; ++i)
      data[i] -= other.data[i];
    return *this;
  }

  matrix& operator*=(const T& scalar)
  {
    size_type n = rows * cols;
    for (size_type i = 0; i < n; ++i)
      data[i] *= scalar;
    return *this;
  }

  matrix& operator/=(const T& scalar)
  {
    if (scalar == T{0})
      throw matrix_exception::not_invertible(); // Reuse for div-by-zero
    T inv = T{1} / scalar;
    size_type n = rows * cols;
    for (size_type i = 0; i < n; ++i)
      data[i] *= inv;
    return *this;
  }

  matrix operator-() const
  {
    matrix res(rows, cols);
    size_type n = rows * cols;
    for (size_type i = 0; i < n; ++i)
      res.data[i] = -data[i];
    return res;
  }

  // Functional application for linear transformation concept
  vector<T, Alloc> operator()(const vector<T, Alloc>& v) const
  {
    return *this * v;
  }

  bool operator==(const matrix& other) const
  {
    if (rows != other.rows || cols != other.cols)
      return false;
    size_type n = rows * cols;
    for (size_type i = 0; i < n; ++i)
      if (data[i] != other.data[i])
        return false;
    return true;
  }
};

// Free function operators for matrix arithmetic

export template<typename T, typename Alloc>
matrix<T, Alloc> operator+(const matrix<T, Alloc>& a, const matrix<T, Alloc>& b)
{
  matrix<T, Alloc> res(a);
  res += b;
  return res;
}

export template<typename T, typename Alloc>
matrix<T, Alloc> operator-(const matrix<T, Alloc>& a, const matrix<T, Alloc>& b)
{
  matrix<T, Alloc> res(a);
  res -= b;
  return res;
}

export template<typename T, typename Alloc>
matrix<T, Alloc> operator*(const T& scalar, const matrix<T, Alloc>& m)
{
  matrix<T, Alloc> res(m);
  res *= scalar;
  return res;
}

export template<typename T, typename Alloc>
matrix<T, Alloc> operator*(const matrix<T, Alloc>& m, const T& scalar)
{
  return scalar * m;
}

export template<typename T, typename Alloc>
matrix<T, Alloc> operator/(const matrix<T, Alloc>& m, const T& scalar)
{
  matrix<T, Alloc> res(m);
  res /= scalar;
  return res;
}

export template<typename T, typename Alloc>
matrix<T, Alloc> operator*(const matrix<T, Alloc>& a, const matrix<T, Alloc>& b)
{
  if (a.get_cols() != b.get_rows())
    throw matrix_exception::dim_mismatch();
  matrix<T, Alloc> res(a.get_rows(), b.get_cols());
  for (std::size_t i = 0; i < a.get_rows(); ++i)
  {
    for (std::size_t k = 0; k < a.get_cols(); ++k)
    {
      for (std::size_t j = 0; j < b.get_cols(); ++j)
      {
        res[i, j] += a[i, k] * b[k, j];
      }
    }
  }
  return res;
}

export template<typename T, typename Alloc>
vector<T, Alloc> operator*(const matrix<T, Alloc>& m, const vector<T, Alloc>& v)
{
  if (m.get_cols() != v.size())
    throw matrix_exception::dim_mismatch();
  vector<T, Alloc> res(m.get_rows());
  for (std::size_t i = 0; i < m.get_rows(); ++i)
  {
    for (std::size_t j = 0; j < m.get_cols(); ++j)
    {
      res[i] += m[i, j] * v[j];
    }
  }
  return res;
}

// Verify matrix satisfies concepts
  static_assert(lam::linalg::concepts::experimental::matrix_c_weak<matrix<double>, double>);
  static_assert(lam::linalg::concepts::experimental::linear_transformation_c_weak<matrix<double>, vector<double>, vector<double>>);

} // namespace lam::linalg

export template<typename T, typename Alloc>
struct std::formatter<lam::linalg::matrix<T, Alloc>>
{
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  auto format(const lam::linalg::matrix<T, Alloc>& m, std::format_context& ctx) const
  {
    auto out = ctx.out();
    std::format_to(out, "[\n");
    for (std::size_t i = 0; i < m.get_rows(); ++i)
    {
      std::format_to(out, "  (");
      for (std::size_t j = 0; j < m.get_cols(); ++j)
      {
        if (j > 0)
          std::format_to(out, ", ");
        std::format_to(out, "{}", m[i, j]);
      }
      std::format_to(out, ")\n");
    }
    return std::format_to(out, "]");
  }
};
