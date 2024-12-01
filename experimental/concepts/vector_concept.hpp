
#include <concepts>
#include <memory>
#include <numeric>
#include <type_traits>

namespace linalg
{

namespace concepts
{

template<typename T>
concept arithmetic = std::is_arithmetic<T>::value;

template<typename G>
concept group_element_c_weak = requires(G g, G h, G k)
{
  { -g } -> std::same_as<G>;
  { g + h } -> std::same_as<decltype(h + g)>;
  { g - h } -> std::same_as<decltype(h - g)>;
  { (g + h) + k } -> std::same_as<decltype(g + (h + k))>;
};

template<typename V>
concept vector_c = requires(V u, V v)
{
  { u + v } -> std::same_as<decltype(v + u)>;
  { u - v } -> std::same_as<decltype(v - u)>;
};

} // end namespace concept
} // end namespace linalg
