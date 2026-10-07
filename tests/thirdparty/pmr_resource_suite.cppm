/*
 *  pmr_resource_suite.cppm - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  Module linalg.test.pmr_suite: one suite, run against any
 *  std::pmr::memory_resource. Shared by the Boost and TBB third-party tests so
 *  both ask the same questions of their resources, and so a std resource can
 *  serve as the control in either file. Imports std + lam.linearalgebra ONLY --
 *  no Boost, no TBB; each test TU includes its own third-party headers.
 */

export module linalg.test.pmr_suite;

import std;
import lam.linearalgebra;

export {
  // Returns the number of failed checks, printing each. The caller owns its own failure count.
  int run_resource_suite(std::string_view name, std::pmr::memory_resource* mr)
  {
    int failures = 0;
    auto check = [&](bool ok, std::string_view what) {
      if (!ok)
      {
        std::println("FAILED [{}]: {}", name, what);
        ++failures;
      }
    };

    lam::linalg::pmr::vector<double> v({3.0, 4.0, 0.0}, mr);
    check(v.size() == 3, "vector size");
    check(v.norm() == 5.0, "vector norm");

    lam::linalg::pmr::vector<double> a({1.0, 2.0}, mr);
    lam::linalg::pmr::vector<double> b({3.0, 4.0}, mr);
    auto sum = a + b;
    check(sum[1] == 6.0, "operator+ value");
    check(sum.get_allocator() == a.get_allocator(), "operator+ result changed resource");

    lam::linalg::pmr::matrix<double> m(8, 8, mr);
    for (std::size_t i = 0; i < 8; ++i)
      m[i, i] = 2.0;
    auto product = m * m;
    check(product[3, 3] == 4.0, "operator* value");
    check(product.get_allocator() == m.get_allocator(), "operator* result changed resource");

    // Create and destroy repeatedly: an arena only grows, a pool recycles, and both must behave.
    for (int i = 0; i < 16; ++i)
    {
      lam::linalg::pmr::vector<double> scratch(64, mr);
      scratch[0] = static_cast<double>(i);
      check(scratch[0] == static_cast<double>(i), "vector over repeated create/destroy");
    }

    if (failures == 0)
      std::println("{} passed", name);
    return failures;
  }
}
