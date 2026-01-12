/*
 *  vectorspace_exceptions.cppm - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  vectorspace_exceptions is a c++ module
 */

export module lam.linearalgebra:vectorspace.exceptions;

import std;

namespace lam::linalg
{

export struct vector_exception : public std::exception
{
  struct non_pos : public std::exception
  {
    const char* what() const noexcept override { return "The dimension of the vector must be at least 1.\n"; }
  };
  struct div_by_zero : public std::exception
  {
    const char* what() const noexcept override { return "Division by zero is undefined!\n"; }
  };
  struct cross_undef : public std::exception
  {
    const char* what() const noexcept override { return "Cross product only defined in 3D.\n"; }
  };
  struct out_of_bounds : public std::exception
  {
    const char* what() const noexcept override { return "Index out of bounds.\n"; }
  };
  struct null_construction : public std::exception
  {
    const char* what() const noexcept override { return "Construction with nullptr not allowed.\n"; }
  };
  struct plus_equals_unequal_dim : public std::exception
  {
    const char* what() const noexcept override { return "Dimension mismatch for operator+=.\n"; }
  };
};

} // namespace lam::linalg
