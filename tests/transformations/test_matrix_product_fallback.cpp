#include <cassert>
#include <iostream>

import lam.linearalgebra;

int main()
{
  using namespace lam::linalg;

  // Test row_major * row_major fallback
  matrix<int> A_rr(2, 3);
  matrix<int> B_rr(3, 2);
  A_rr[0, 0] = 1;
  A_rr[0, 1] = 2;
  A_rr[0, 2] = 3;
  A_rr[1, 0] = 4;
  A_rr[1, 1] = 5;
  A_rr[1, 2] = 6;
  B_rr[0, 0] = 10;
  B_rr[0, 1] = 11;
  B_rr[1, 0] = 20;
  B_rr[1, 1] = 21;
  B_rr[2, 0] = 30;
  B_rr[2, 1] = 31;

  auto C1 = A_rr * B_rr;
  assert((C1[0, 0] == 140));
  assert((C1[0, 1] == 146));
  assert((C1[1, 0] == 320));
  assert((C1[1, 1] == 335));

  // Test col_major * col_major fallback
  matrix<int, std::allocator<int>, storage_layout::col_major> A_cc(2, 3);
  matrix<int, std::allocator<int>, storage_layout::col_major> B_cc(3, 2);
  A_cc[0, 0] = 1;
  A_cc[0, 1] = 2;
  A_cc[0, 2] = 3;
  A_cc[1, 0] = 4;
  A_cc[1, 1] = 5;
  A_cc[1, 2] = 6;
  B_cc[0, 0] = 10;
  B_cc[0, 1] = 11;
  B_cc[1, 0] = 20;
  B_cc[1, 1] = 21;
  B_cc[2, 0] = 30;
  B_cc[2, 1] = 31;

  auto C2 = A_cc * B_cc;
  assert((C2[0, 0] == 140));
  assert((C2[0, 1] == 146));
  assert((C2[1, 0] == 320));
  assert((C2[1, 1] == 335));

  // Test col_major * row_major fallback (uses generic tiled cache approach)
  matrix<int, std::allocator<int>, storage_layout::col_major> A_cr(2, 3);
  matrix<int, std::allocator<int>, storage_layout::row_major> B_cr(3, 2);
  A_cr[0, 0] = 1;
  A_cr[0, 1] = 2;
  A_cr[0, 2] = 3;
  A_cr[1, 0] = 4;
  A_cr[1, 1] = 5;
  A_cr[1, 2] = 6;
  B_cr[0, 0] = 10;
  B_cr[0, 1] = 11;
  B_cr[1, 0] = 20;
  B_cr[1, 1] = 21;
  B_cr[2, 0] = 30;
  B_cr[2, 1] = 31;

  auto C3 = A_cr * B_cr;
  assert((C3[0, 0] == 140));
  assert((C3[0, 1] == 146));
  assert((C3[1, 0] == 320));
  assert((C3[1, 1] == 335));

  // Test row_major * col_major fallback
  matrix<int, std::allocator<int>, storage_layout::row_major> A_rc(2, 3);
  matrix<int, std::allocator<int>, storage_layout::col_major> B_rc(3, 2);
  A_rc[0, 0] = 1;
  A_rc[0, 1] = 2;
  A_rc[0, 2] = 3;
  A_rc[1, 0] = 4;
  A_rc[1, 1] = 5;
  A_rc[1, 2] = 6;
  B_rc[0, 0] = 10;
  B_rc[0, 1] = 11;
  B_rc[1, 0] = 20;
  B_rc[1, 1] = 21;
  B_rc[2, 0] = 30;
  B_rc[2, 1] = 31;

  auto C4 = A_rc * B_rc;
  assert((C4[0, 0] == 140));
  assert((C4[0, 1] == 146));
  assert((C4[1, 0] == 320));
  assert((C4[1, 1] == 335));

  std::cout << "All fallback matrix multiplications passed." << std::endl;
  return 0;
}
