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
  struct unequal_allocators : public std::exception
  {
    const char* what() const noexcept override { return "Swap of unequal allocators that do not propagate."; }
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

inline std::optional<int> checked_lwork(const double query, const int info)
{
  if (info != 0)
    return std::nullopt;
  if (std::isnan(query) || query < 1.0 || query > static_cast<double>(std::numeric_limits<int>::max()))
    return std::nullopt;
  return static_cast<int>(query);
}
} // namespace detail

export template<typename T, storage_layout Layout = storage_layout::row_major, typename Alloc = std::allocator<T>>
  requires lam::concepts::experimental::ring_element_c_weak<T>
class matrix;

export template<typename T, storage_layout Layout, typename Alloc>
  requires lam::concepts::experimental::ring_element_c_weak<T>
class matrix
{
  template<typename U, storage_layout L, typename A>
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

  constexpr explicit matrix(const Alloc& alloc) noexcept : m_rows{0}, m_cols{0}, m_storage(alloc) {}

  constexpr matrix(size_type r, size_type c)
    requires std::default_initializable<Alloc>
    : matrix(r, c, Alloc())
  {}

  constexpr matrix(size_type r, size_type c, const Alloc& alloc) : m_rows{r}, m_cols{c}, m_storage(r * c, alloc)
  {
    if (r == 0 || c == 0)
      throw matrix_exception::non_pos();
  }

  constexpr explicit matrix(size_type n)
    requires std::default_initializable<Alloc>
    : matrix(n, Alloc())
  {}

  constexpr matrix(size_type n, const Alloc& alloc) : matrix(n, n, alloc)
  {
    for (size_type i = 0; i < n; ++i)
      (*this)[i, i] = T{1};
  }

  constexpr matrix(const std::vector<std::vector<T>>& mtrx)
    requires std::default_initializable<Alloc>
    : matrix(mtrx, Alloc())
  {}

  constexpr matrix(const std::vector<std::vector<T>>& mtrx, const Alloc& alloc)
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

  constexpr matrix(matrix&& other) noexcept
    : m_rows{std::exchange(other.m_rows, 0)}, m_cols{std::exchange(other.m_cols, 0)},
      m_storage{std::move(other.m_storage)}, label{std::move(other.label)}
  {}

  constexpr matrix& operator=(matrix&& other) noexcept(std::is_nothrow_move_assignable_v<vector<T, Alloc>>)
  {
    if (this == &other)
      return *this;
    m_storage = std::move(other.m_storage);
    m_rows = other.m_rows;
    m_cols = other.m_cols;
    if (other.m_storage.size() == 0)
    {
      other.m_rows = 0;
      other.m_cols = 0;
    }
    label = std::move(other.label);
    return *this;
  }

  matrix(const matrix& other) = default;
  matrix& operator=(const matrix& other) = default;

  // Conversion constructor from matrix with different layout
  template<storage_layout OtherLayout>
  constexpr matrix(const matrix<T, OtherLayout, Alloc>& other)
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

  constexpr auto transpose() const
  {
    // Return opposite layout: RowMajor -> ColMajor, ColMajor -> RowMajor
    constexpr storage_layout NewLayout =
      (Layout == storage_layout::row_major) ? storage_layout::col_major : storage_layout::row_major;

    matrix<T, NewLayout, Alloc> res(m_cols, m_rows, get_allocator());

    // row_major(r, c) and col_major(c, r) hold the same linear data, so the transpose is a flat copy.
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

  template<typename P>
  static constexpr auto row_of(P base, size_type r, size_type c, size_type i)
  {
    if constexpr (Layout == storage_layout::row_major)
      return std::span{base + i * c, c};
    else
      return std::span{base + i, (c - 1) * r + 1} | std::views::stride(r);
  }

  template<typename P>
  static constexpr auto col_of(P base, size_type r, size_type c, size_type j)
  {
    if constexpr (Layout == storage_layout::col_major)
      return std::span{base + j * r, r};
    else
      return std::span{base + j, (r - 1) * c + 1} | std::views::stride(c);
  }

  void row(size_type) const&& = delete;
  constexpr auto row(size_type i) &
  {
    if (i >= m_rows)
      throw matrix_exception::out_of_bounds();
    return row_of(m_storage.begin(), m_rows, m_cols, i);
  }

  constexpr auto row(size_type i) const&
  {
    if (i >= m_rows)
      throw matrix_exception::out_of_bounds();
    return row_of(m_storage.begin(), m_rows, m_cols, i);
  }

  void rows_range() const&& = delete;
  constexpr auto rows_range() &
  {
    return std::views::iota(size_type{0}, m_rows) |
           std::views::transform(
             [base = m_storage.begin(), r = m_rows, c = m_cols](size_type i) { return row_of(base, r, c, i); });
  }

  constexpr auto rows_range() const&
  {
    return std::views::iota(size_type{0}, m_rows) |
           std::views::transform(
             [base = m_storage.begin(), r = m_rows, c = m_cols](size_type i) { return row_of(base, r, c, i); });
  }

  void col(size_type) const&& = delete;
  constexpr auto col(size_type j) &
  {
    if (j >= m_cols)
      throw matrix_exception::out_of_bounds();
    return col_of(m_storage.begin(), m_rows, m_cols, j);
  }

  constexpr auto col(size_type j) const&
  {
    if (j >= m_cols)
      throw matrix_exception::out_of_bounds();
    return col_of(m_storage.begin(), m_rows, m_cols, j);
  }

  void cols_range() const&& = delete;
  constexpr auto cols_range() &
  {
    return std::views::iota(size_type{0}, m_cols) |
           std::views::transform(
             [base = m_storage.begin(), r = m_rows, c = m_cols](size_type j) { return col_of(base, r, c, j); });
  }

  constexpr auto cols_range() const&
  {
    return std::views::iota(size_type{0}, m_cols) |
           std::views::transform(
             [base = m_storage.begin(), r = m_rows, c = m_cols](size_type j) { return col_of(base, r, c, j); });
  }

  template<typename A2>
  constexpr matrix& operator+=(const matrix<T, Layout, A2>& other)
  {
    if (m_rows != other.m_rows || m_cols != other.m_cols)
      throw matrix_exception::dim_mismatch();
    m_storage += other.m_storage;
    return *this;
  }

  template<typename A2>
  constexpr matrix& operator-=(const matrix<T, Layout, A2>& other)
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

  // Throws when the allocators differ and do not propagate; swap never copies elements.
  constexpr void swap(matrix& other) noexcept(std::is_nothrow_swappable_v<vector<T, Alloc>>)
  {
    if constexpr (!std::allocator_traits<Alloc>::propagate_on_container_swap::value)
      if (!(get_allocator() == other.get_allocator()))
        throw matrix_exception::unequal_allocators();
    m_storage.swap(other.m_storage);
    std::swap(m_rows, other.m_rows);
    std::swap(m_cols, other.m_cols);
    label.swap(other.label);
  }
  friend constexpr void swap(matrix& a, matrix& b) noexcept(noexcept(a.swap(b))) { a.swap(b); }

  static constexpr matrix zero()
    requires std::default_initializable<Alloc>
  { return matrix{}; }

  constexpr bool operator==(const matrix& other) const
  {
    if (m_rows != other.m_rows || m_cols != other.m_cols)
      return false;
    return m_storage == other.m_storage;
  }

  constexpr void set_label(std::string_view l) { label = l; }
  constexpr std::optional<std::string> get_label() const { return label; }

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

export template<typename T, storage_layout LayoutA, storage_layout LayoutB, typename A1, typename A2>
constexpr matrix<T, LayoutA, A1> operator+(const matrix<T, LayoutA, A1>& a, const matrix<T, LayoutB, A2>& b)
{
  matrix<T, LayoutA, A1> res(a);
  if constexpr (LayoutA == LayoutB)
    res += b;
  else
  {
    if (a.rows() != b.rows() || a.cols() != b.cols())
      throw matrix_exception::dim_mismatch();
    for (std::size_t i = 0; i < a.rows(); ++i)
      for (std::size_t j = 0; j < a.cols(); ++j)
        res[i, j] += b[i, j];
  }
  return res;
}

export template<typename T, storage_layout LayoutA, storage_layout LayoutB, typename A1, typename A2>
constexpr matrix<T, LayoutA, A1> operator-(const matrix<T, LayoutA, A1>& a, const matrix<T, LayoutB, A2>& b)
{
  matrix<T, LayoutA, A1> res(a);
  if constexpr (LayoutA == LayoutB)
    res -= b;
  else
  {
    if (a.rows() != b.rows() || a.cols() != b.cols())
      throw matrix_exception::dim_mismatch();
    for (std::size_t i = 0; i < a.rows(); ++i)
      for (std::size_t j = 0; j < a.cols(); ++j)
        res[i, j] -= b[i, j];
  }
  return res;
}

export template<typename T, storage_layout Layout, typename Alloc>
constexpr matrix<T, Layout, Alloc> operator*(const T& scalar, const matrix<T, Layout, Alloc>& m)
{
  matrix<T, Layout, Alloc> res(m);
  res *= scalar;
  return res;
}

export template<typename T, storage_layout Layout, typename Alloc>
constexpr matrix<T, Layout, Alloc> operator*(const matrix<T, Layout, Alloc>& m, const T& scalar)
{ return scalar * m; }

export template<typename T, storage_layout Layout, typename Alloc>
constexpr matrix<T, Layout, Alloc> operator/(const matrix<T, Layout, Alloc>& m, const T& scalar)
{
  matrix<T, Layout, Alloc> res(m);
  res /= scalar;
  return res;
}


export template<typename T, storage_layout Layout, typename Alloc>
constexpr matrix<T, Layout, Alloc> zero_like(const matrix<T, Layout, Alloc>& m)
{ return matrix<T, Layout, Alloc>(m.rows(), m.cols(), m.get_allocator()); }

export template<typename T, storage_layout Layout, typename Alloc>
constexpr matrix<T, Layout, Alloc> identity_like(const matrix<T, Layout, Alloc>& m)
{
  if (m.rows() != m.cols())
    throw matrix_exception::dim_mismatch();
  return matrix<T, Layout, Alloc>(m.rows(), m.get_allocator());
}

export template<typename T, storage_layout Layout, typename Alloc>
constexpr vector<T, Alloc> operator*(const matrix<T, Layout, Alloc>& m, const vector<T, Alloc>& v)
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

export namespace pmr
{
template<typename T, storage_layout Layout = storage_layout::row_major>
using matrix = lam::linalg::matrix<T, Layout, std::pmr::polymorphic_allocator<T>>;
} // namespace pmr

} // namespace lam::linalg

export template<typename T, lam::linalg::storage_layout Layout, typename Alloc>
struct std::formatter<lam::linalg::matrix<T, Layout, Alloc>>
{
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  auto format(const lam::linalg::matrix<T, Layout, Alloc>& m, std::format_context& ctx) const
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
