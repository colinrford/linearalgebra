/*
 *  matrix.cppm - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  matrix is a c++ module
 */

export module lam.linearalgebra:matrix;
import std;
import lam.concepts;
import :config;
import :vectorspace;

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
    const char* what() const noexcept override { return "Index out of matrix bounds."; }
  };
  struct not_invertible : public std::exception
  {
    const char* what() const noexcept override { return "Matrix is not invertible."; }
  };
};

export enum class storage_layout
{
  row_major,
  col_major
};

namespace detail
{
// Allocator for a T-valued result derived from M: M's allocator rebound to T, else std::allocator<T>.
template<typename M, typename T>
struct result_alloc
{
  using type = std::allocator<T>;
};

template<typename M, typename T>
  requires requires { typename M::allocator_type; }
struct result_alloc<M, T>
{
  using type = typename std::allocator_traits<typename M::allocator_type>::template rebind_alloc<T>;
};

template<typename M, typename T = typename M::scalar_type>
using result_alloc_t = typename result_alloc<M, T>::type;

template<typename M, typename T = typename M::scalar_type>
constexpr result_alloc_t<M, T> result_allocator(const M& m)
{
  if constexpr (requires { m.get_allocator(); })
    return result_alloc_t<M, T>(m.get_allocator());
  else
    return result_alloc_t<M, T>();
}
} // namespace detail

export template<typename T, typename Alloc = std::allocator<T>, storage_layout Layout = storage_layout::row_major>
  requires lam::concepts::experimental::ring_element_c_weak<T>
class matrix;

export template<typename T, typename Alloc, storage_layout Layout>
  requires lam::concepts::experimental::ring_element_c_weak<T>
class matrix
{
  template<typename U, typename A, storage_layout L>
    requires lam::concepts::experimental::ring_element_c_weak<U>
  friend class matrix;

public:
  using value_type = T;
  using scalar_type = T;
  using allocator_type = Alloc;
  using size_type = std::size_t;
  static constexpr storage_layout layout = Layout;

private:
  size_type m_rows;
  size_type m_cols;
  vector<T, Alloc> m_storage;
  std::optional<std::string> label;

  constexpr size_type index(size_type i, size_type j) const noexcept
  {
    if constexpr (Layout == storage_layout::row_major)
      return i * m_cols + j;
    else
      return j * m_rows + i;
  }

public:
  constexpr matrix() : m_rows{0}, m_cols{0}, m_storage{} {}

  constexpr matrix(size_type r, size_type c, const Alloc& alloc = Alloc())
    : m_rows{r}, m_cols{c}, m_storage(r * c, alloc)
  {
    if (r == 0 || c == 0)
      throw matrix_exception::non_pos();
  }

  constexpr explicit matrix(size_type n, const Alloc& alloc = Alloc()) : matrix(n, n, alloc)
  {
    for (size_type i = 0; i < n; ++i)
      (*this)[i, i] = T{1};
  }

  matrix(const std::vector<std::vector<T>>& mtrx, const Alloc& alloc = Alloc())
    : matrix(mtrx.size(), mtrx.empty() ? 0 : mtrx[0].size(), alloc)
  {
    for (size_type i = 0; i < m_rows; ++i)
    {
      for (size_type j = 0; j < m_cols; ++j)
      {
        (*this)[i, j] = mtrx[i][j];
      }
    }
  }

  matrix(matrix&&) noexcept = default;
  matrix& operator=(matrix&&) noexcept = default;

  matrix(const matrix& other) = default;
  matrix& operator=(const matrix& other) = default;

  // Conversion constructor from matrix with different layout
  template<storage_layout OtherLayout>
  constexpr matrix(const matrix<T, Alloc, OtherLayout>& other)
    : m_rows{other.rows()}, m_cols{other.cols()}, m_storage(other.rows() * other.cols(), other.get_allocator())
  {
    // Shuffle copy
    for (size_type i = 0; i < m_rows; ++i)
    {
      for (size_type j = 0; j < m_cols; ++j)
      {
        (*this)[i, j] = other[i, j];
      }
    }
  }
  // C++23 Multidimensional []
  constexpr T& operator[](size_type i, size_type j)
  {
    if (i >= m_rows || j >= m_cols)
      throw matrix_exception::out_of_bounds();
    return m_storage[index(i, j)];
  }

  constexpr const T& operator[](size_type i, size_type j) const
  {
    if (i >= m_rows || j >= m_cols)
      throw matrix_exception::out_of_bounds();
    return m_storage[index(i, j)];
  }

  constexpr size_type rows() const noexcept { return m_rows; }
  constexpr size_type cols() const noexcept { return m_cols; }

  auto transpose() const
  {
    // Return opposite layout: RowMajor -> ColMajor, ColMajor -> RowMajor
    constexpr storage_layout NewLayout =
      (Layout == storage_layout::row_major) ? storage_layout::col_major : storage_layout::row_major;

    matrix<T, Alloc, NewLayout> res(m_cols, m_rows, get_allocator());

    // row_major(r, c) and col_major(c, r) hold the same linear data, so the transpose is a flat copy.
    size_type n = m_rows * m_cols;
    std::copy(m_storage.begin(), m_storage.end(), res.begin());

    return res;
  }

  // ========== Iterators and Accessors ==========

  // Element iteration (row-major order)
  constexpr T* begin() noexcept { return m_storage.begin(); }
  constexpr T* end() noexcept { return m_storage.end(); }
  constexpr const T* begin() const noexcept { return m_storage.begin(); }
  constexpr const T* end() const noexcept { return m_storage.end(); }
  constexpr const T* cbegin() const noexcept { return m_storage.cbegin(); }
  constexpr const T* cend() const noexcept { return m_storage.cend(); }

  // Single row access — returns contiguous span
  auto row(size_type i)
  {
    if (i >= m_rows)
      throw matrix_exception::out_of_bounds();

    if constexpr (Layout == storage_layout::row_major)
    {
      return std::span<T>{m_storage.begin() + i * m_cols, m_cols};
    }
    else
    {
      return std::views::iota(size_type{0}, m_cols) |
             std::views::transform([this, i](size_type j) -> T& { return (*this)[i, j]; });
    }
  }

  auto row(size_type i) const
  {
    if (i >= m_rows)
      throw matrix_exception::out_of_bounds();

    if constexpr (Layout == storage_layout::row_major)
    {
      return std::span<const T>{m_storage.begin() + i * m_cols, m_cols};
    }
    else
    {
      return std::views::iota(size_type{0}, m_cols) |
             std::views::transform([this, i](size_type j) -> const T& { return (*this)[i, j]; });
    }
  }

  // All rows — returns range of spans
  auto rows_range()
  {
    return std::views::iota(size_type{0}, m_rows) | std::views::transform([this](size_type i) { return this->row(i); });
  }

  auto rows_range() const
  {
    return std::views::iota(size_type{0}, m_rows) | std::views::transform([this](size_type i) { return this->row(i); });
  }

  // Single column access — strided view (non-contiguous)
  auto col(size_type j)
  {
    if (j >= m_cols)
      throw matrix_exception::out_of_bounds();

    if constexpr (Layout == storage_layout::col_major)
    {
      return std::span<T>{m_storage.begin() + j * m_rows, m_rows};
    }
    else
    {
      return std::views::iota(size_type{0}, m_rows) |
             std::views::transform([this, j](size_type i) -> T& { return (*this)[i, j]; });
    }
  }

  auto col(size_type j) const
  {
    if (j >= m_cols)
      throw matrix_exception::out_of_bounds();

    if constexpr (Layout == storage_layout::col_major)
    {
      return std::span<const T>{m_storage.begin() + j * m_rows, m_rows};
    }
    else
    {
      return std::views::iota(size_type{0}, m_rows) |
             std::views::transform([this, j](size_type i) -> const T& { return (*this)[i, j]; });
    }
  }

  // All columns — returns range of strided ranges
  auto cols_range()
  {
    return std::views::iota(size_type{0}, m_cols) | std::views::transform([this](size_type j) { return this->col(j); });
  }

  auto cols_range() const
  {
    return std::views::iota(size_type{0}, m_cols) | std::views::transform([this](size_type j) { return this->col(j); });
  }

  constexpr matrix& operator+=(const matrix& other)
  {
    if (m_rows != other.m_rows || m_cols != other.m_cols)
      throw matrix_exception::dim_mismatch();
    m_storage += other.m_storage;
    return *this;
  }

  constexpr matrix& operator-=(const matrix& other)
  {
    if (m_rows != other.m_rows || m_cols != other.m_cols)
      throw matrix_exception::dim_mismatch();
    m_storage -= other.m_storage;
    return *this;
  }

  constexpr matrix& operator*=(const T& scalar)
  {
    m_storage *= scalar;
    return *this;
  }

  constexpr matrix& operator/=(const T& scalar)
  {
    m_storage /= scalar;
    return *this;
  }

  constexpr matrix operator-() const
  {
    matrix res(*this);
    for (auto& x : res.m_storage)
      x = -x;
    return res;
  }

  // Functional application for linear transformation concept
  constexpr vector<T, Alloc> operator()(const vector<T, Alloc>& v) const { return *this * v; }

  constexpr allocator_type get_allocator() const noexcept { return m_storage.get_allocator(); }

  bool operator==(const matrix& other) const
  {
    if (m_rows != other.m_rows || m_cols != other.m_cols)
      return false;
    return m_storage == other.m_storage;
  }

  void set_label(std::string_view l) { label = l; }
  std::optional<std::string> get_label() const { return label; }

  // ========== Query Methods ==========

  constexpr bool is_square() const noexcept { return m_rows == m_cols; }

  constexpr T trace() const
  {
    if (!is_square())
      throw matrix_exception::dim_mismatch();
    T sum = T{0};
    for (size_type i = 0; i < m_rows; ++i)
      sum += (*this)[i, i];
    return sum;
  }

  constexpr vector<T, Alloc> diagonal() const
  {
    size_type diag_size = std::min(m_rows, m_cols);
    vector<T, Alloc> diag(diag_size, get_allocator());
    for (size_type i = 0; i < diag_size; ++i)
      diag[i] = (*this)[i, i];
    return diag;
  }
};

// Free function operators for matrix arithmetic

export template<typename T, typename Alloc, storage_layout Layout>
constexpr matrix<T, Alloc, Layout> operator+(const matrix<T, Alloc, Layout>& a, const matrix<T, Alloc, Layout>& b)
{
  matrix<T, Alloc, Layout> res(a);
  res += b;
  return res;
}

export template<typename T, typename Alloc, storage_layout Layout>
constexpr matrix<T, Alloc, Layout> operator-(const matrix<T, Alloc, Layout>& a, const matrix<T, Alloc, Layout>& b)
{
  matrix<T, Alloc, Layout> res(a);
  res -= b;
  return res;
}

export template<typename T, typename Alloc, storage_layout Layout>
constexpr matrix<T, Alloc, Layout> operator*(const T& scalar, const matrix<T, Alloc, Layout>& m)
{
  matrix<T, Alloc, Layout> res(m);
  res *= scalar;
  return res;
}

export template<typename T, typename Alloc, storage_layout Layout>
constexpr matrix<T, Alloc, Layout> operator*(const matrix<T, Alloc, Layout>& m, const T& scalar)
{ return scalar * m; }

export template<typename T, typename Alloc, storage_layout Layout>
constexpr matrix<T, Alloc, Layout> operator/(const matrix<T, Alloc, Layout>& m, const T& scalar)
{
  matrix<T, Alloc, Layout> res(m);
  res /= scalar;
  return res;
}


export template<typename T, typename Alloc, storage_layout Layout>
constexpr vector<T, Alloc> operator*(const matrix<T, Alloc, Layout>& m, const vector<T, Alloc>& v)
{
  if (m.cols() != v.size())
    throw matrix_exception::dim_mismatch();
  vector<T, Alloc> res(m.rows(), m.get_allocator());
  for (std::size_t i = 0; i < m.rows(); ++i)
  {
    for (std::size_t j = 0; j < m.cols(); ++j)
    {
      res[i] += m[i, j] * v[j];
    }
  }
  return res;
}

// Verify matrix satisfies concepts
static_assert(lam::concepts::experimental::matrix_c_weak<matrix<double>, double>);
static_assert(
  lam::concepts::experimental::linear_transformation_c_weak<matrix<double>, vector<double>, vector<double>>);

} // namespace lam::linalg

export template<typename T, typename Alloc, lam::linalg::storage_layout Layout>
struct std::formatter<lam::linalg::matrix<T, Alloc, Layout>>
{
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  auto format(const lam::linalg::matrix<T, Alloc, Layout>& m, std::format_context& ctx) const
  {
    auto out = ctx.out();
    std::format_to(out, "[\n");
    for (std::size_t i = 0; i < m.rows(); ++i)
    {
      std::format_to(out, "  (");
      for (std::size_t j = 0; j < m.cols(); ++j)
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
