/*
 *  linalg.cppm - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  linalg is a c++ module
 */

export module lam.linearalgebra;
export import :vectorspace;
export import lam.concepts;
export import :matrix;
export import :matrix.product;
export import :matrix.decomposition;
export import :matrix.eigenvalue;
export import :matrix.svd;

export import :config;

export namespace lam
{
using lam::linalg::matrix;
using lam::linalg::vector;
} // namespace lam
