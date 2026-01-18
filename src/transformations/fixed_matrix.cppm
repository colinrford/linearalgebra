/*
 *  fixed_matrix.cppm - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  fixed_matrix: A stack-allocated, fixed-size matrix type.
 *  Column-major storage for GPU compatibility.
 *  Designed for constexpr contexts and graphics/physics applications.
 *
 *  NOTE: This partition is NOT automatically exported.
 *  To enable, add to linalg.cppm:
 *    export import :fixed_matrix;
 */

module;

import std;

export module lam.linearalgebra:fixed_matrix;

import lam.concepts;
import :vectorspace.fixed_vector;

namespace lam::linalg
{

// ============================================================================
// fixed_matrix<T, Rows, Cols> - Stack-allocated, fixed-size matrix
// Column-major storage (matches OpenGL/Vulkan/WebGPU conventions)
// ============================================================================
export template<typename T, std::size_t Rows, std::size_t Cols>
  requires lam::concepts::experimental::ring_element_c_weak<T>
struct fixed_matrix
{
  using scalar_type = T;
  using size_type = std::size_t;
  using value_type = T;

  // Column-major: data[col * Rows + row]
  std::array<T, Rows * Cols> data{};

  // ========== Static Info ==========

  static constexpr size_type static_rows = Rows;
  static constexpr size_type static_cols = Cols;

  // ========== Constructors ==========

  constexpr fixed_matrix() = default;

  // From flat array (column-major order)
  constexpr explicit fixed_matrix(const std::array<T, Rows * Cols>& arr) noexcept
    : data{arr}
  {}

  // From initializer (column-major order)
  template<typename... Args>
    requires (sizeof...(Args) == Rows * Cols) && (std::convertible_to<Args, T> && ...)
  constexpr fixed_matrix(Args... args) noexcept
    : data{static_cast<T>(args)...}
  {}

  // ========== Dimensions ==========

  constexpr size_type rows() const noexcept { return Rows; }
  constexpr size_type cols() const noexcept { return Cols; }
  constexpr bool is_square() const noexcept { return Rows == Cols; }

  // ========== Element Access (column-major) ==========

  constexpr T& operator[](size_type row, size_type col) noexcept
  {
    return data[col * Rows + row];
  }

  constexpr const T& operator[](size_type row, size_type col) const noexcept
  {
    return data[col * Rows + row];
  }

  constexpr T& at(size_type row, size_type col)
  {
    if (row >= Rows || col >= Cols)
      throw std::out_of_range("fixed_matrix index out of bounds");
    return data[col * Rows + row];
  }

  constexpr const T& at(size_type row, size_type col) const
  {
    if (row >= Rows || col >= Cols)
      throw std::out_of_range("fixed_matrix index out of bounds");
    return data[col * Rows + row];
  }

  // ========== Column Access ==========

  constexpr fixed_vector<T, Rows> col(size_type c) const noexcept
  {
    fixed_vector<T, Rows> result;
    for (size_type r = 0; r < Rows; ++r)
      result[r] = (*this)[r, c];
    return result;
  }

  constexpr void set_col(size_type c, const fixed_vector<T, Rows>& v) noexcept
  {
    for (size_type r = 0; r < Rows; ++r)
      (*this)[r, c] = v[r];
  }

  // ========== Row Access ==========

  constexpr fixed_vector<T, Cols> row(size_type r) const noexcept
  {
    fixed_vector<T, Cols> result;
    for (size_type c = 0; c < Cols; ++c)
      result[c] = (*this)[r, c];
    return result;
  }

  // ========== Static Factories ==========

  static constexpr fixed_matrix identity() noexcept
    requires (Rows == Cols)
  {
    fixed_matrix m{};
    for (size_type i = 0; i < Rows; ++i)
      m[i, i] = T{1};
    return m;
  }

  static constexpr fixed_matrix zero() noexcept
  {
    return fixed_matrix{};
  }

  // ========== Arithmetic ==========

  constexpr fixed_matrix operator-() const noexcept
  {
    fixed_matrix result;
    for (size_type i = 0; i < Rows * Cols; ++i)
      result.data[i] = -data[i];
    return result;
  }

  constexpr fixed_matrix& operator+=(const fixed_matrix& other) noexcept
  {
    for (size_type i = 0; i < Rows * Cols; ++i)
      data[i] += other.data[i];
    return *this;
  }

  constexpr fixed_matrix& operator-=(const fixed_matrix& other) noexcept
  {
    for (size_type i = 0; i < Rows * Cols; ++i)
      data[i] -= other.data[i];
    return *this;
  }

  constexpr fixed_matrix& operator*=(const T& scalar) noexcept
  {
    for (size_type i = 0; i < Rows * Cols; ++i)
      data[i] *= scalar;
    return *this;
  }

  constexpr fixed_matrix& operator/=(const T& scalar) noexcept
  {
    for (size_type i = 0; i < Rows * Cols; ++i)
      data[i] /= scalar;
    return *this;
  }

  // ========== Matrix-Vector Multiplication ==========

  constexpr fixed_vector<T, Rows> operator*(const fixed_vector<T, Cols>& v) const noexcept
  {
    fixed_vector<T, Rows> result{};
    for (size_type r = 0; r < Rows; ++r)
    {
      for (size_type c = 0; c < Cols; ++c)
      {
        result[r] += (*this)[r, c] * v[c];
      }
    }
    return result;
  }

  // ========== Transpose ==========

  constexpr fixed_matrix<T, Cols, Rows> transpose() const noexcept
  {
    fixed_matrix<T, Cols, Rows> result;
    for (size_type r = 0; r < Rows; ++r)
    {
      for (size_type c = 0; c < Cols; ++c)
      {
        result[c, r] = (*this)[r, c];
      }
    }
    return result;
  }

  // ========== Trace (square matrices only) ==========

  constexpr T trace() const noexcept
    requires (Rows == Cols)
  {
    T sum = T{0};
    for (size_type i = 0; i < Rows; ++i)
      sum += (*this)[i, i];
    return sum;
  }

  // ========== Comparison ==========

  constexpr bool operator==(const fixed_matrix&) const noexcept = default;

  // ========== Raw Data Access (for GPU upload) ==========

  constexpr const T* raw() const noexcept { return data.data(); }
  constexpr T* raw() noexcept { return data.data(); }

  // ========== Iterators (column-major order) ==========

  constexpr T* begin() noexcept { return data.data(); }
  constexpr T* end() noexcept { return data.data() + Rows * Cols; }
  constexpr const T* begin() const noexcept { return data.data(); }
  constexpr const T* end() const noexcept { return data.data() + Rows * Cols; }
};

// ============================================================================
// Matrix-Matrix Multiplication
// ============================================================================

export template<typename T, std::size_t M, std::size_t N, std::size_t P>
constexpr fixed_matrix<T, M, P> operator*(const fixed_matrix<T, M, N>& a,
                                           const fixed_matrix<T, N, P>& b) noexcept
{
  fixed_matrix<T, M, P> result{};
  for (std::size_t r = 0; r < M; ++r)
  {
    for (std::size_t c = 0; c < P; ++c)
    {
      for (std::size_t k = 0; k < N; ++k)
      {
        result[r, c] += a[r, k] * b[k, c];
      }
    }
  }
  return result;
}

// ============================================================================
// Free Function Operators
// ============================================================================

export template<typename T, std::size_t Rows, std::size_t Cols>
constexpr fixed_matrix<T, Rows, Cols> operator+(const fixed_matrix<T, Rows, Cols>& a,
                                                 const fixed_matrix<T, Rows, Cols>& b) noexcept
{
  fixed_matrix<T, Rows, Cols> result = a;
  result += b;
  return result;
}

export template<typename T, std::size_t Rows, std::size_t Cols>
constexpr fixed_matrix<T, Rows, Cols> operator-(const fixed_matrix<T, Rows, Cols>& a,
                                                 const fixed_matrix<T, Rows, Cols>& b) noexcept
{
  fixed_matrix<T, Rows, Cols> result = a;
  result -= b;
  return result;
}

export template<typename T, std::size_t Rows, std::size_t Cols>
constexpr fixed_matrix<T, Rows, Cols> operator*(const T& scalar,
                                                 const fixed_matrix<T, Rows, Cols>& m) noexcept
{
  fixed_matrix<T, Rows, Cols> result = m;
  result *= scalar;
  return result;
}

export template<typename T, std::size_t Rows, std::size_t Cols>
constexpr fixed_matrix<T, Rows, Cols> operator*(const fixed_matrix<T, Rows, Cols>& m,
                                                 const T& scalar) noexcept
{
  return scalar * m;
}

export template<typename T, std::size_t Rows, std::size_t Cols>
constexpr fixed_matrix<T, Rows, Cols> operator/(const fixed_matrix<T, Rows, Cols>& m,
                                                 const T& scalar) noexcept
{
  fixed_matrix<T, Rows, Cols> result = m;
  result /= scalar;
  return result;
}

// ============================================================================
// 4x4 Matrix Utilities (Common for Graphics)
// ============================================================================

export template<typename T>
constexpr fixed_matrix<T, 4, 4> translate(const fixed_vector<T, 3>& t) noexcept
{
  auto m = fixed_matrix<T, 4, 4>::identity();
  m[0, 3] = t.x();
  m[1, 3] = t.y();
  m[2, 3] = t.z();
  return m;
}

export template<typename T>
constexpr fixed_matrix<T, 4, 4> scale(const fixed_vector<T, 3>& s) noexcept
{
  fixed_matrix<T, 4, 4> m{};
  m[0, 0] = s.x();
  m[1, 1] = s.y();
  m[2, 2] = s.z();
  m[3, 3] = T{1};
  return m;
}

export template<typename T>
constexpr fixed_matrix<T, 4, 4> scale(T s) noexcept
{
  return scale(fixed_vector<T, 3>{s, s, s});
}

export template<typename T>
constexpr fixed_matrix<T, 4, 4> rotate_x(T angle) noexcept
{
  T c = std::cos(angle);
  T s = std::sin(angle);
  auto m = fixed_matrix<T, 4, 4>::identity();
  m[1, 1] = c;  m[1, 2] = -s;
  m[2, 1] = s;  m[2, 2] = c;
  return m;
}

export template<typename T>
constexpr fixed_matrix<T, 4, 4> rotate_y(T angle) noexcept
{
  T c = std::cos(angle);
  T s = std::sin(angle);
  auto m = fixed_matrix<T, 4, 4>::identity();
  m[0, 0] = c;   m[0, 2] = s;
  m[2, 0] = -s;  m[2, 2] = c;
  return m;
}

export template<typename T>
constexpr fixed_matrix<T, 4, 4> rotate_z(T angle) noexcept
{
  T c = std::cos(angle);
  T s = std::sin(angle);
  auto m = fixed_matrix<T, 4, 4>::identity();
  m[0, 0] = c;  m[0, 1] = -s;
  m[1, 0] = s;  m[1, 1] = c;
  return m;
}

// Rotation around arbitrary axis (Rodrigues' formula)
export template<typename T>
constexpr fixed_matrix<T, 4, 4> rotate(T angle, const fixed_vector<T, 3>& axis) noexcept
{
  // Normalize axis
  T len = std::sqrt(axis.x()*axis.x() + axis.y()*axis.y() + axis.z()*axis.z());
  if (len == T{0}) return fixed_matrix<T, 4, 4>::identity();

  T x = axis.x() / len;
  T y = axis.y() / len;
  T z = axis.z() / len;

  T c = std::cos(angle);
  T s = std::sin(angle);
  T t = T{1} - c;

  fixed_matrix<T, 4, 4> m{};
  m[0, 0] = t*x*x + c;     m[0, 1] = t*x*y - s*z;  m[0, 2] = t*x*z + s*y;  m[0, 3] = T{0};
  m[1, 0] = t*x*y + s*z;   m[1, 1] = t*y*y + c;    m[1, 2] = t*y*z - s*x;  m[1, 3] = T{0};
  m[2, 0] = t*x*z - s*y;   m[2, 1] = t*y*z + s*x;  m[2, 2] = t*z*z + c;    m[2, 3] = T{0};
  m[3, 0] = T{0};          m[3, 1] = T{0};         m[3, 2] = T{0};         m[3, 3] = T{1};
  return m;
}

// Look-at matrix (view matrix)
export template<typename T>
constexpr fixed_matrix<T, 4, 4> look_at(const fixed_vector<T, 3>& eye,
                                         const fixed_vector<T, 3>& center,
                                         const fixed_vector<T, 3>& up) noexcept
{
  // Forward (camera looks down -Z in view space)
  fixed_vector<T, 3> f{
    center.x() - eye.x(),
    center.y() - eye.y(),
    center.z() - eye.z()
  };
  T f_len = std::sqrt(f.x()*f.x() + f.y()*f.y() + f.z()*f.z());
  if (f_len > T{0}) { f.x() /= f_len; f.y() /= f_len; f.z() /= f_len; }

  // Right = f x up
  fixed_vector<T, 3> r{
    f.y()*up.z() - f.z()*up.y(),
    f.z()*up.x() - f.x()*up.z(),
    f.x()*up.y() - f.y()*up.x()
  };
  T r_len = std::sqrt(r.x()*r.x() + r.y()*r.y() + r.z()*r.z());
  if (r_len > T{0}) { r.x() /= r_len; r.y() /= r_len; r.z() /= r_len; }

  // Up = r x f
  fixed_vector<T, 3> u{
    r.y()*f.z() - r.z()*f.y(),
    r.z()*f.x() - r.x()*f.z(),
    r.x()*f.y() - r.y()*f.x()
  };

  // Dot products for translation
  T dot_r_eye = r.x()*eye.x() + r.y()*eye.y() + r.z()*eye.z();
  T dot_u_eye = u.x()*eye.x() + u.y()*eye.y() + u.z()*eye.z();
  T dot_f_eye = f.x()*eye.x() + f.y()*eye.y() + f.z()*eye.z();

  fixed_matrix<T, 4, 4> m{};
  m[0, 0] = r.x();  m[0, 1] = r.y();  m[0, 2] = r.z();   m[0, 3] = -dot_r_eye;
  m[1, 0] = u.x();  m[1, 1] = u.y();  m[1, 2] = u.z();   m[1, 3] = -dot_u_eye;
  m[2, 0] = -f.x(); m[2, 1] = -f.y(); m[2, 2] = -f.z();  m[2, 3] = dot_f_eye;
  m[3, 0] = T{0};   m[3, 1] = T{0};   m[3, 2] = T{0};    m[3, 3] = T{1};
  return m;
}

// Perspective projection matrix
export template<typename T>
constexpr fixed_matrix<T, 4, 4> perspective(T fov_y_radians, T aspect,
                                             T near, T far) noexcept
{
  T tan_half = std::tan(fov_y_radians / T{2});
  T range = far - near;

  fixed_matrix<T, 4, 4> m{};
  m[0, 0] = T{1} / (aspect * tan_half);
  m[1, 1] = T{1} / tan_half;
  m[2, 2] = -(far + near) / range;
  m[2, 3] = -T{2} * far * near / range;
  m[3, 2] = -T{1};
  // m[3, 3] = 0 (default)
  return m;
}

// Orthographic projection matrix
export template<typename T>
constexpr fixed_matrix<T, 4, 4> ortho(T left, T right, T bottom, T top,
                                       T near, T far) noexcept
{
  fixed_matrix<T, 4, 4> m{};
  m[0, 0] = T{2} / (right - left);
  m[1, 1] = T{2} / (top - bottom);
  m[2, 2] = -T{2} / (far - near);
  m[0, 3] = -(right + left) / (right - left);
  m[1, 3] = -(top + bottom) / (top - bottom);
  m[2, 3] = -(far + near) / (far - near);
  m[3, 3] = T{1};
  return m;
}

// ============================================================================
// Common Type Aliases
// ============================================================================

export using mat2f = fixed_matrix<float, 2, 2>;
export using mat3f = fixed_matrix<float, 3, 3>;
export using mat4f = fixed_matrix<float, 4, 4>;

export using mat2d = fixed_matrix<double, 2, 2>;
export using mat3d = fixed_matrix<double, 3, 3>;
export using mat4d = fixed_matrix<double, 4, 4>;

export using mat2x3f = fixed_matrix<float, 2, 3>;
export using mat3x2f = fixed_matrix<float, 3, 2>;
export using mat2x4f = fixed_matrix<float, 2, 4>;
export using mat4x2f = fixed_matrix<float, 4, 2>;
export using mat3x4f = fixed_matrix<float, 3, 4>;
export using mat4x3f = fixed_matrix<float, 4, 3>;

} // namespace lam::linalg

// ============================================================================
// std::formatter specialization
// ============================================================================

export template<typename T, std::size_t Rows, std::size_t Cols>
struct std::formatter<lam::linalg::fixed_matrix<T, Rows, Cols>>
{
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  auto format(const lam::linalg::fixed_matrix<T, Rows, Cols>& m,
              std::format_context& ctx) const
  {
    auto out = ctx.out();
    std::format_to(out, "[\n");
    for (std::size_t r = 0; r < Rows; ++r)
    {
      std::format_to(out, "  [");
      for (std::size_t c = 0; c < Cols; ++c)
      {
        if (c > 0) std::format_to(out, ", ");
        std::format_to(out, "{:8.4f}", m[r, c]);
      }
      std::format_to(out, "]\n");
    }
    return std::format_to(out, "]");
  }
};
