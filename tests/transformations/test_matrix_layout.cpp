import std;
import lam.linearalgebra;
import lam.concepts;

using namespace lam::linalg;

int main()
{
  constexpr auto row_layout = storage_layout::row_major;
  constexpr auto col_layout = storage_layout::col_major;

  // 1. Storage Verification
  // Create a 2x3 matrix
  // [ 1 2 3 ]
  // [ 4 5 6 ]

  matrix<int, std::allocator<int>, row_layout> m_row(2, 3);
  matrix<int, std::allocator<int>, col_layout> m_col(2, 3);

  auto filler = [](auto& m) {
    m[0, 0] = 1;
    m[0, 1] = 2;
    m[0, 2] = 3;
    m[1, 0] = 4;
    m[1, 1] = 5;
    m[1, 2] = 6;
  };

  filler(m_row);
  filler(m_col);

  // Check logical equality
  if (m_row == m_row)
    std::println("m_row equality checked."); // Just to define operator usage

  // Row Major: m[1, 0] is element idx 3 (1*3 + 0) -> Data: 1, 2, 3, 4, 5, 6
  // Col Major: m[1, 0] is element idx 1 (0*2 + 1) -> Wait: index = j*rows + i
  // Element (1,0): i=1, j=0. Index = 0*2 + 1 = 1.
  // Element (0,1): i=0, j=1. Index = 1*2 + 0 = 2.

  // Verify underlying data
  // m_row data: [1, 2, 3, 4, 5, 6]
  // m_col data:
  // (0,0)=1 -> idx 0
  // (1,0)=4 -> idx 1
  // (0,1)=2 -> idx 2
  // (1,1)=5 -> idx 3
  // (0,2)=3 -> idx 4
  // (1,2)=6 -> idx 5
  // So m_col data should be: [1, 4, 2, 5, 3, 6]

  auto check_data = [](auto& m, std::vector<int> expected) {
    int idx = 0;
    for (auto val : m)
    { // Iterates underlying storage
      if (val != expected[idx])
      {
        std::println("Mismatch at index {}: got {}, expected {}", idx, val, expected[idx]);
        std::exit(1);
      }
      idx++;
    }
  };

  std::println("Verifying Row Major Storage...");
  check_data(m_row, {1, 2, 3, 4, 5, 6});
  std::println("PASSED");

  std::println("Verifying Column Major Storage...");
  check_data(m_col, {1, 4, 2, 5, 3, 6});
  std::println("PASSED");

  // 2. Algorithm Verification
  // Check if generic algorithms accept col_major matrix
  matrix<double, std::allocator<double>, col_layout> m_double(2, 2);
  m_double[0, 0] = 4.0;
  m_double[0, 1] = 3.0;
  m_double[1, 0] = 6.0;
  m_double[1, 1] = 3.0;
  // Det = 12 - 18 = -6

  auto det = crout_lu_det(m_double);
  if (std::abs(det - (-6.0)) < 1e-9)
  {
    std::println("Algorthm Generic Check: crout_lu_det passed for Col Major. Det: {}", det);
  }
  else
  {
    std::println("Algorthm Generic Check: FAILED. Det: {}", det);
    return 1;
  }

  // 3. Row/Col View Types
  // Row major: row() is span, col() is view
  static_assert(std::is_same_v<decltype(m_row.row(0)), std::span<int>>);
  static_assert(!std::is_same_v<decltype(m_row.col(0)), std::span<int>>);

  // Col major: row() is view, col() is span
  static_assert(!std::is_same_v<decltype(m_col.row(0)), std::span<int>>);
  static_assert(std::is_same_v<decltype(m_col.col(0)), std::span<int>>);

  std::println("Type System Checks: PASSED");

  // 4. Mixed Layout Multiplication
  // m_mix_row (RowMajor) * m_mix_col (ColMajor)
  matrix<int, std::allocator<int>, row_layout> m_mix_row(2, 2);
  m_mix_row[0, 0] = 1;
  m_mix_row[0, 1] = 2;
  m_mix_row[1, 0] = 3;
  m_mix_row[1, 1] = 4;

  matrix<int, std::allocator<int>, col_layout> m_mix_col(2, 2);
  m_mix_col[0, 0] = 1;
  m_mix_col[0, 1] = 0;
  m_mix_col[1, 0] = 0;
  m_mix_col[1, 1] = 1;

  // Identity multiplication should be unchanged
  auto m_res = m_mix_row * m_mix_col;

  if (m_res[0, 0] == 1 && m_res[0, 1] == 2 && m_res[1, 0] == 3 && m_res[1, 1] == 4)
  {
    std::println("Mixed Layout Mul (Row * Col): PASSED");
  }
  else
  {
    std::println("Mixed Layout Mul (Row * Col): FAILED");
    return 1;
  }

  // 5. Transpose Optimization
  // RowMajor transpose -> ColMajor (O(N) copy)
  auto tr = m_mix_row.transpose();
  static_assert(decltype(tr)::layout == storage_layout::col_major, "Transpose should switch layout");

  // Logical check
  // m_mix_row:
  // 1 2
  // 3 4
  // tr (ColMajor):
  // 1 3 (col 0: 1, 3)
  // 2 4 (col 1: 2, 4)
  if (tr[0, 0] == 1 && tr[0, 1] == 3 && tr[1, 0] == 2 && tr[1, 1] == 4)
  {
    std::println("Transpose layout switch values: PASSED");
  }
  else
  {
    std::println("Transpose layout switch values: FAILED");
    return 1;
  }

  // Explicit conversion back
  matrix<int, std::allocator<int>, row_layout> back_to_row = tr;
  if (back_to_row == matrix<int, std::allocator<int>, row_layout>(2, 2))
  {
    // Actually verify values
    if (back_to_row[0, 1] == 3) // A^T (row major)
      std::println("Conversion back to RowMajor: PASSED");
  }

  return 0;
}
