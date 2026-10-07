/*
 *  test_accessor_lifetime.cpp - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  test_accessor_lifetime is a c++ module
 */

import std;
import lam.linearalgebra;

using lam::linalg::matrix;
using lam::linalg::storage_layout;

int failures = 0;

void check(bool ok, std::string_view what)
{
  if (!ok)
  {
    std::println("FAILED: {}", what);
    ++failures;
  }
}

template<storage_layout Layout>
matrix<double, Layout> make()
{ return matrix<double, Layout>({{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}}); }

template<typename M>
concept row_callable = requires(M&& m) { std::forward<M>(m).row(std::size_t{0}); };
template<typename M>
concept col_callable = requires(M&& m) { std::forward<M>(m).col(std::size_t{0}); };
template<typename M>
concept rows_range_callable = requires(M&& m) { std::forward<M>(m).rows_range(); };
template<typename M>
concept cols_range_callable = requires(M&& m) { std::forward<M>(m).cols_range(); };

template<storage_layout Layout>
void lifetime_suite(std::string_view name)
{
  using mat = matrix<double, Layout>;

  static_assert(row_callable<mat&> && row_callable<const mat&>);
  static_assert(col_callable<mat&> && col_callable<const mat&>);
  static_assert(rows_range_callable<mat&> && rows_range_callable<const mat&>);
  static_assert(cols_range_callable<mat&> && cols_range_callable<const mat&>);

  static_assert(!row_callable<mat> && !row_callable<const mat>);
  static_assert(!col_callable<mat> && !col_callable<const mat>);
  static_assert(!rows_range_callable<mat> && !rows_range_callable<const mat>);
  static_assert(!cols_range_callable<mat> && !cols_range_callable<const mat>);

  const auto m = make<Layout>();

  check(std::ranges::equal(m.row(1), std::vector<double>{4.0, 5.0, 6.0}), std::format("{}: row values", name));
  check(std::ranges::equal(m.col(2), std::vector<double>{3.0, 6.0}), std::format("{}: col values", name));
  check(std::ranges::size(m.rows_range()) == 2, std::format("{}: rows_range size", name));
  check(std::ranges::size(m.cols_range()) == 3, std::format("{}: cols_range size", name));

  double sum = 0.0;
  for (auto r : m.rows_range())
    for (const double& x : r)
      sum += x;
  check(sum == 21.0, std::format("{}: rows_range values", name));

  sum = 0.0;
  for (auto c : m.cols_range())
    for (const double& x : c)
      sum += x;
  check(sum == 21.0, std::format("{}: cols_range values", name));
}

// Both layouts must agree here: an accessor follows the storage it was built from, not the matrix object.
template<storage_layout Layout>
void moved_from_suite(std::string_view name)
{
  auto m = make<Layout>();
  auto r = m.row(0);
  auto moved = std::move(m);
  check(std::ranges::equal(r, std::vector<double>{1.0, 2.0, 3.0}), std::format("{}: row survives a move", name));
}

consteval bool constexpr_accessors()
{
  const matrix<double, storage_layout::col_major> m({{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}});
  double sum = 0.0;
  for (const double& x : m.row(1))
    sum += x;
  return sum == 15.0 && m.col(2)[1] == 6.0;
}

static_assert(constexpr_accessors());

int main()
{
  lifetime_suite<storage_layout::row_major>("row_major");
  lifetime_suite<storage_layout::col_major>("col_major");

  moved_from_suite<storage_layout::row_major>("row_major");
  moved_from_suite<storage_layout::col_major>("col_major");

  if (failures != 0)
  {
    std::println("{} failure(s)", failures);
    return 1;
  }
  std::println("\n✅ matrix accessor lifetime tests passed");
  return 0;
}
