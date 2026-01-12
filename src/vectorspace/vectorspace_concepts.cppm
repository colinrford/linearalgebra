/*
 *  vectorspace_concepts.cppm - Colin Ford
 *    see github.com/colinrford/linearalgebra for more info
 *    lam.linearalgebra is unlicensed at this time
 *
 *  vectorspace_concepts is a c++ module
 */

export module lam.linearalgebra:vectorspace.concepts;

import std;
import lam.concepts;

namespace lam::linalg
{

namespace concepts::experimental::internals
{
/*
 *  left_module_element_c_weak
 *  Models a left module M over a scalar ring R.
 *  Scalar multiplication: r * m -> M
 */
template<typename M, typename R>
concept left_module_element_c_weak = lam::concepts::experimental::additive_group_element_c_weak<M> and
                                     lam::concepts::experimental::ring_element_c_weak<R> and requires(M m, R r) {
                                       { r * m } -> std::same_as<M>;
                                     };

/*
 *  right_module_element_c_weak
 *  Models a right module M over a scalar ring R.
 *  Scalar multiplication: m * r -> M
 */
template<typename M, typename R>
concept right_module_element_c_weak = lam::concepts::experimental::additive_group_element_c_weak<M> and
                                      lam::concepts::experimental::ring_element_c_weak<R> and requires(M m, R r) {
                                        { m * r } -> std::same_as<M>;
                                      };

/*
 *  module_element_c_weak
 *  Models a module M that supports both left and right scalar multiplication.
 *  (Often effectively a bimodule where the left and right actions are compatible,
 *   or a module over a commutative ring where left and right are the same).
 */
template<typename M, typename R>
concept module_element_c_weak = left_module_element_c_weak<M, R> and right_module_element_c_weak<M, R>;

/*
 *  vectorspace_element_c_weak
 *  Models a vector space V over a scalar field K.
 *  Requirements:
 *  1. It is a module over K (supporting both left/right mult for this library's conventions).
 *  2. K is a field.
 */
template<typename V, typename K>
concept vectorspace_element_c_weak =
  module_element_c_weak<V, K> and lam::concepts::experimental::field_element_c_weak<K>;

/*
 *  vector_c_weak
 *  Models a concrete vector container (like lam::linalg::vector).
 */
template<typename V>
concept vector_c_weak =
  requires {
    typename V::scalar_type;
    typename V::size_type;
  } and module_element_c_weak<V, typename V::scalar_type> and
  requires(V v, V::size_type i) {
    { v.size() } -> std::same_as<typename V::size_type>;
    { v[i] } -> std::convertible_to<typename V::scalar_type>;
  } and
  std::same_as<std::remove_cvref_t<decltype(std::declval<V>()[std::declval<typename V::size_type>()])>,
               typename V::scalar_type>;
} // namespace concepts::experimental::internals

namespace concepts::experimental
{
export template<typename M, typename R>
concept left_module_element_c_weak = internals::left_module_element_c_weak<M, R>;

export template<typename M, typename R>
concept right_module_element_c_weak = internals::right_module_element_c_weak<M, R>;

export template<typename M, typename R>
concept module_element_c_weak = internals::module_element_c_weak<M, R>;

export template<typename V, typename K>
concept vectorspace_element_c_weak = internals::vectorspace_element_c_weak<V, K>;

export template<typename V>
concept vector_c_weak = internals::vector_c_weak<V>;
} // namespace concepts::experimental

} // namespace lam::linalg
