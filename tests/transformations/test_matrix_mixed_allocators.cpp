/*
 *  test_matrix_mixed_allocators.cpp - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  test_matrix_mixed_allocators is a c++ module
 */

import std;
import lam.linearalgebra;

using lam::linalg::matrix;
using lam::linalg::storage_layout;

template<typename T>
struct tagged_alloc
{
  using value_type = T;
  int tag = 0;

  constexpr tagged_alloc() = default;
  constexpr explicit tagged_alloc(int t) : tag{t} {}
  template<typename U>
  constexpr tagged_alloc(const tagged_alloc<U>& o) : tag{o.tag}
  {}

  constexpr T* allocate(std::size_t n) { return std::allocator<T>{}.allocate(n); }
  constexpr void deallocate(T* p, std::size_t n) { std::allocator<T>{}.deallocate(p, n); }

  constexpr bool operator==(const tagged_alloc& o) const { return tag == o.tag; }
};

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
using std_m = matrix<double, Layout, std::allocator<double>>;
template<storage_layout Layout>
using tag_m = matrix<double, Layout, tagged_alloc<double>>;
template<storage_layout Layout>
using pmr_m = matrix<double, Layout, std::pmr::polymorphic_allocator<double>>;

const std::vector<std::vector<double>> lhs_data{{1.0, 2.0}, {3.0, 4.0}};
const std::vector<std::vector<double>> rhs_data{{10.0, 20.0}, {30.0, 40.0}};

template<storage_layout Layout>
void mixed_allocator_suite(std::string_view name)
{
  std::array<std::byte, 8192> buffer{};
  std::pmr::monotonic_buffer_resource pool{buffer.data(), buffer.size()};

  const std_m<Layout> a(lhs_data);
  const tag_m<Layout> b(rhs_data, tagged_alloc<double>{7});
  const pmr_m<Layout> c(rhs_data, &pool);

  static_assert(std::same_as<decltype(a + b), std_m<Layout>>);
  static_assert(std::same_as<decltype(b + a), tag_m<Layout>>);
  static_assert(std::same_as<decltype(c - a), pmr_m<Layout>>);

  {
    auto s = a + b;
    check(s[0, 0] == 11.0 && s[1, 1] == 44.0, std::format("{}: mixed operator+ values", name));
  }
  {
    auto s = b + a;
    check(s[0, 0] == 11.0 && s[1, 1] == 44.0, std::format("{}: mixed operator+ values (reversed)", name));
    check(s.get_allocator().tag == 7, std::format("{}: operator+ must keep the left operand's allocator", name));
  }
  {
    auto d = c - a;
    check(d[0, 0] == 9.0 && d[1, 1] == 36.0, std::format("{}: mixed operator- values", name));
    check(d.get_allocator().resource() == &pool,
          std::format("{}: operator- must keep the left operand's resource", name));
  }
}

template<storage_layout LayoutA, storage_layout LayoutB>
void mixed_product_suite(std::string_view name)
{
  std::array<std::byte, 8192> buffer{};
  std::pmr::monotonic_buffer_resource pool{buffer.data(), buffer.size()};

  const tag_m<LayoutA> a(lhs_data, tagged_alloc<double>{3});
  const pmr_m<LayoutB> b(rhs_data, &pool);

  static_assert(std::same_as<decltype(a * b), matrix<double, storage_layout::row_major, tagged_alloc<double>>>);

  auto p = a * b;
  check(p[0, 0] == 70.0 && p[1, 1] == 220.0, std::format("{}: mixed operator* values", name));
  check(p.get_allocator().tag == 3, std::format("{}: operator* must keep the left operand's allocator", name));
}

const std::vector<std::vector<double>> wide_lhs{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};
const std::vector<std::vector<double>> wide_rhs{{10.0, 20.0, 30.0}, {40.0, 50.0, 60.0}};

template<typename M>
bool equals_rows(const M& m, const std::vector<std::vector<double>>& expected)
{
  for (std::size_t i = 0; i < expected.size(); ++i)
    for (std::size_t j = 0; j < expected[i].size(); ++j)
      if (m[i, j] != expected[i][j])
        return false;
  return true;
}

template<storage_layout LayoutA, storage_layout LayoutB>
void mixed_layout_suite(std::string_view name)
{
  std::array<std::byte, 8192> buffer{};
  std::pmr::monotonic_buffer_resource pool{buffer.data(), buffer.size()};

  const std_m<LayoutA> a(wide_lhs);
  const std_m<LayoutB> b(wide_rhs);
  const tag_m<LayoutA> t(wide_lhs, tagged_alloc<double>{5});
  const pmr_m<LayoutB> p(wide_rhs, &pool);

  static_assert(std::same_as<decltype(a + b), std_m<LayoutA>>);
  static_assert(std::same_as<decltype(b - a), std_m<LayoutB>>);
  static_assert(std::same_as<decltype(t + p), tag_m<LayoutA>>);
  static_assert(std::same_as<decltype(p - t), pmr_m<LayoutB>>);

  check(equals_rows(a + b, {{11.0, 22.0, 33.0}, {44.0, 55.0, 66.0}}), std::format("{}: operator+ values", name));
  check(equals_rows(b - a, {{9.0, 18.0, 27.0}, {36.0, 45.0, 54.0}}), std::format("{}: operator- values", name));
  {
    auto s = t + p;
    check(equals_rows(s, {{11.0, 22.0, 33.0}, {44.0, 55.0, 66.0}}), std::format("{}: mixed operator+ values", name));
    check(s.get_allocator().tag == 5, std::format("{}: operator+ must keep the left operand's allocator", name));
  }
  {
    auto d = p - t;
    check(equals_rows(d, {{9.0, 18.0, 27.0}, {36.0, 45.0, 54.0}}), std::format("{}: mixed operator- values", name));
    check(d.get_allocator().resource() == &pool,
          std::format("{}: operator- must keep the left operand's resource", name));
  }

  const std_m<LayoutB> tall({{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}});
  bool threw = false;
  try
  {
    static_cast<void>(a + tall);
  }
  catch (const lam::linalg::matrix_exception::dim_mismatch&)
  {
    threw = true;
  }
  check(threw, std::format("{}: 2x3 + 3x2 must throw dim_mismatch", name));
}

template<typename X, typename Y>
concept addable = requires(const X& x, const Y& y) { x + y; };

static_assert(addable<std_m<storage_layout::row_major>, tag_m<storage_layout::row_major>>);
static_assert(
  !addable<std_m<storage_layout::row_major>, matrix<float, storage_layout::row_major, std::allocator<float>>>);

consteval bool constexpr_mixed()
{
  const std::vector<std::vector<double>> l{{1.0, 2.0}, {3.0, 4.0}};
  const std::vector<std::vector<double>> r{{10.0, 20.0}, {30.0, 40.0}};

  const std_m<storage_layout::row_major> a(l);
  const tag_m<storage_layout::row_major> b(r, tagged_alloc<double>{1});

  auto s = a + b;
  if (s[0, 0] != 11.0 || s[1, 1] != 44.0)
    return false;

  auto p = b * a;
  if (p.get_allocator().tag != 1 || p[0, 0] != 70.0)
    return false;

  const matrix<double, storage_layout::col_major, std::allocator<double>> c(r);
  auto m = a + c;
  return m[0, 1] == 22.0 && m[1, 0] == 33.0;
}

static_assert(constexpr_mixed());

int main()
{
  mixed_allocator_suite<storage_layout::row_major>("row_major");
  mixed_allocator_suite<storage_layout::col_major>("col_major");

  mixed_product_suite<storage_layout::row_major, storage_layout::row_major>("row x row");
  mixed_product_suite<storage_layout::row_major, storage_layout::col_major>("row x col");
  mixed_product_suite<storage_layout::col_major, storage_layout::col_major>("col x col");

  mixed_layout_suite<storage_layout::row_major, storage_layout::col_major>("row + col");
  mixed_layout_suite<storage_layout::col_major, storage_layout::row_major>("col + row");
  mixed_layout_suite<storage_layout::row_major, storage_layout::row_major>("row + row");

  if (failures != 0)
  {
    std::println("{} failure(s)", failures);
    return 1;
  }
  std::println("\n✅ matrix mixed-allocator tests passed");
  return 0;
}
