
module;

/*
 *  if stuck with gcc 14
 */
import <concepts>;
import <numeric>;
import <type_traits>;

/*
 *  clang 19 and beyond, gcc 15 and beyond
 */
//import std;

export module la.concepts;

export namespace la
{

namespace concepts
{

/*
 *  all concepts will be declared with _c to distinguish them from types, etc
 */

template<typename T>
concept arithmetic = std::is_arithmetic<T>::value;

/*
 *  operator overloaded weak syntactical 'requirements' modeling group
 *  elements
 */
template<typename G>
concept group_element_c_weak = requires(G g, G h, G k)
{
  { -g };// -> std::same_as<decltype(g)>;
  { g + h } -> std::same_as<decltype(h + g)>;
  { g - h } -> std::same_as<decltype(h - g)>;
  { (g + h) + k } -> std::same_as<decltype(g + (h + k))>;
};

/*
 *  operator overloaded weak syntactical 'requirements' modeling field
 *  elements
 */
template<typename K>
concept field_element_c_weak = group_element_c_weak<K> and
requires(K a, K b, K c)
{
  { 1 / a } -> std::same_as<K>;
  { a * b } -> std::same_as<decltype(b * a)>;
  { (a * b) * c } -> std::same_as<decltype(a * (b * c))>;
};

template<typename V, typename K>
concept vector_c = group_element_c_weak<V> and field_element_c_weak<K>;

} // end namespace concept
} // end namespace la
