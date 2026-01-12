/*
 *  transformations_concepts.cppm - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  transformations_concepts is a c++ module
 */

export module lam.linearalgebra:transformations.concepts;

import std;
import lam.concepts;
import :vectorspace.concepts; // For vectorspace_element_c_weak

namespace lam::linalg::concepts::experimental::internals
{

/**
 * Models a linear transformation T: V -> W.
 * A linear transformation is a map between two vector spaces that preserves
 * vector addition and scalar multiplication.
 */
template<typename F, typename V, typename W>
concept linear_transformation_c_weak =
  lam::linalg::concepts::experimental::vectorspace_element_c_weak<V, typename V::scalar_type> &&
  lam::linalg::concepts::experimental::vectorspace_element_c_weak<W, typename W::scalar_type> && requires(F f, V v) {
    { f(v) } -> std::convertible_to<W>; // Functional application returns vector in W
  };

/**
 * Models a matrix (finitary linear transformation).
 * Supports element access, row/column dimensions, and row/column iterators.
 */
template<typename M, typename Scalar>
concept matrix_c_weak =
  lam::linalg::concepts::experimental::vectorspace_element_c_weak<M, Scalar> // Matrix itself is a vector space
  && requires(M m, typename M::size_type i) {
       typename M::size_type;

       // Dimensions
       // Dimensions
       { m.rows() } -> std::integral;
       { m.cols() } -> std::integral;

       // Element access
       { m[i, i] } -> std::convertible_to<Scalar>;

       // Row/Column accessors (return iterable views)
       { m.row(i) };
       { m.col(i) };
       // { m.rows_range() }; // Optional?
       // { m.cols_range() }; // Optional?
     };

/**
 * Models a representation of a group G on a vector space V.
 * A representation is a homomorphism rho: G -> GL(V).
 * For each g in G, rho(g) is an invertible linear transformation (automorphism) of V.
 */
template<typename Rho, typename G, typename V>
concept representation_c_weak = requires(Rho rho, G g, V v) {
  // rho(g) must be a linear transformation V -> V
  { rho(g) } -> linear_transformation_c_weak<V, V>;

  // Applying the transformation to a vector
  { rho(g)(v) } -> std::same_as<V>;
};

/**
 * Models a square matrix (rows == cols).
 * Structural check: type must have is_square() method.
 */
template<typename M, typename Scalar>
concept square_matrix_c_weak = matrix_c_weak<M, Scalar> && requires(M m) {
  { m.is_square() } -> std::same_as<bool>;
};

/**
 * Models an invertible matrix type.
 * Structural check: type must have det() and inverse() methods.
 * Runtime singularity check still required.
 */
template<typename M, typename Scalar>
concept invertible_c_weak = square_matrix_c_weak<M, Scalar> && requires(M m) {
  { m.det() } -> std::convertible_to<Scalar>;
  { m.inverse() } -> matrix_c_weak<Scalar>;
};

} // namespace lam::linalg::concepts::experimental::internals

namespace lam::linalg::concepts::experimental
{
export template<typename F, typename V, typename W>
concept linear_transformation_c_weak = internals::linear_transformation_c_weak<F, V, W>;

export template<typename M, typename Scalar>
concept matrix_c_weak = internals::matrix_c_weak<M, Scalar>;

export template<typename M, typename Scalar>
concept square_matrix_c_weak = internals::square_matrix_c_weak<M, Scalar>;

export template<typename M, typename Scalar>
concept invertible_c_weak = internals::invertible_c_weak<M, Scalar>;

export template<typename Rho, typename G, typename V>
concept representation_c_weak = internals::representation_c_weak<Rho, G, V>;
} // namespace lam::linalg::concepts::experimental
