
#include <chrono>
#include <list>
#include <numeric>
#include <ratio>
#include <sstream>
#include <string>
#include <thread>
//#include <tbb/tbb.h>
#include "vector.cpp"

template <typename T>
using std_vector = std::vector<T>;
using la_vector = linalg::vector;

//std::tuple<std_vector<la_vector>,
//           std::chrono::duration<double, std::nano>>
template <typename T = la_vector>
std::tuple<std_vector<T>, std::chrono::duration<double, std::nano>>
time_this_process(auto process, auto... args)
{
  const auto t_0 = std::chrono::steady_clock::now();
  std_vector<T> results = process(std::forward<decltype(args)>(args)...);
  const auto t_1 = std::chrono::steady_clock::now();
  const std::chrono::duration<double, std::nano> duration =
        std::chrono::duration_cast<std::chrono::nanoseconds>(t_1 - t_0);
  return std::forward_as_tuple(std::move(results), duration);
}

std_vector<la_vector> construct_n_vectors_of_dim(std::size_t n, std::size_t d)
{
  std_vector<la_vector> vectors;
  for (std::size_t i = 0; i != n; ++i)
  { vectors.push_back(la_vector(d)); }
  return vectors;
}

std_vector<la_vector> list_construct_n_vectors(std::size_t n)
{
  std_vector<la_vector> vectors;
  for (std::size_t i = 0; i != n; ++i)
  { vectors.push_back(la_vector{1.,2.,3.,4.}); }
  return vectors;
}

std_vector<la_vector> fill_zero_vectors(std_vector<la_vector> vectors)
{
  double val = 1.;
  std_vector<la_vector> filled_vectors;
  for (std::size_t i = 0; i != vectors.size(); ++i)
  {
    std::fill(vectors[i].begin(), vectors[i].end(), val);
    val = -val;
  }
  return std::forward<decltype(vectors)>(vectors);
}

auto compute_dot_products(std_vector<la_vector>* left_vectors,
                          std_vector<la_vector>* right_vectors)
{
  std_vector<double> dot_products;
  for (std::size_t i = 0; i != left_vectors->size(); ++i)
  { dot_products.push_back((*left_vectors)[i].dot((*right_vectors)[i])); }
  return dot_products;
}

auto compute_better_dot_products(std_vector<la_vector>* left_vectors,
                                 std_vector<la_vector>* right_vectors)
{
  std_vector<double> dot_products;
  for (std::size_t i = 0; i != left_vectors->size(); ++i)
  {
    dot_products.push_back(
                  (*left_vectors)[i].better_dot((*right_vectors)[i])
                );
  }
  return dot_products;
}

/*auto compute_even_better_dot_products(std_vector<la_vector>* left_vectors,
                                      std_vector<la_vector>* right_vectors)
{
  std_vector<double> dot_products;
  for (int i = 0; i < left_vectors->size(); i++)
  {
    dot_products.push_back(
                  (*left_vectors)[i].even_better_dot((*right_vectors)[i])
                );
  }
  return dot_products;
} */

std_vector<std_vector<double>>
construct_n_std_vectors_of_dim(std::size_t n, std::size_t d)
{
  std_vector<std_vector<double>> vectors;
  for (std::size_t i = 0; i < n; i++)
  { vectors.push_back(std_vector<double>(d)); }
  return vectors;
}

std_vector<std_vector<double>> list_construct_n_std_vectors(std::size_t n)
{
  std_vector<std_vector<double>> vectors;
  for (std::size_t i = 0; i < n; i++)
  { vectors.push_back(std_vector<double>{1.,2.,3.,4.}); }
  return vectors;
}

std_vector<std_vector<double>>
fill_zero_std_vectors(std_vector<std_vector<double>> vectors)
{
  double val = 1.;
  std_vector<std_vector<double>> filled_vectors;
  for (std::size_t i = 0; i < vectors.size(); i++)
  {
    std::fill(vectors[i].begin(), vectors[i].end(), val);
    val = -val;
  }
  return std::forward<decltype(vectors)>(vectors);
}

auto compute_std_dot_products(std_vector<std_vector<double>>* left_vectors,
                              std_vector<std_vector<double>>* right_vectors)
{
  std_vector<double> dot_products;
  for (std::size_t i = 0; i < left_vectors->size(); i++)
  {
    dot_products.push_back(std::inner_product((*left_vectors)[i].begin(),
                                              (*left_vectors)[i].end(),
                                              (*right_vectors)[i].begin(),
                                              0.)
                          );
  }
  return dot_products;
}

auto
compute_better_std_dot_products(std_vector<std_vector<double>>* left_vectors,
                                std_vector<std_vector<double>>* right_vectors)
{
  std_vector<double> dot_products;
  for (std::size_t i = 0; i < left_vectors->size(); i++)
  {
    dot_products.push_back(std::transform_reduce((*left_vectors)[i].begin(),
                                                 (*left_vectors)[i].end(),
                                                 (*right_vectors)[i].begin(),
                                                 0.)
                          );
  }
  return dot_products;
}

void print_process_duration(const auto duration,
                            const int number_of_objects,
                            const std::string object_name,
                            const std::string process_name)
{
  std::cout << "Time of process "
            << process_name
            << " for "
            << number_of_objects
            << " "
            << object_name
            << "s took about\t"
            << duration
            << std::endl
            << "\ti.e., " << duration / number_of_objects
            << " per object(s)"
            << std::endl;
}

int main()
{
  std::thread tony;
  const std::string la_vector_name = "la_vector";
  const std::string std_vector_name = "std_vector";
  std::cout << "BEGIN" << std::endl;
  for (int i = 1; i < 100000000; i = i * 10)
  {
    auto [vectors,
          construction_time] = time_this_process(construct_n_vectors_of_dim,
                                                 i, 4);
    print_process_duration(construction_time, i, la_vector_name,
                           "construct_n_vectors_of_dim");
    auto [more_vectors,
          list_construction_time] = time_this_process(list_construct_n_vectors,
                                                      i);
    print_process_duration(list_construction_time, i, la_vector_name,
                           "list_construct_n_vectors");
    auto [filled_vectors,
          fill_time] = time_this_process(fill_zero_vectors,
                                         std::move(vectors));
    print_process_duration(fill_time, i, la_vector_name,
                           "fill_zero_vectors");

    auto [dot_products,
          dot_time] = time_this_process<double>(compute_dot_products,
                                        &filled_vectors, &more_vectors);

    print_process_duration(dot_time, i, la_vector_name,
                           "compute_dot_products");

    auto [better_dot_products,
          better_dot_time] = time_this_process<double>
                             (compute_better_dot_products,
                             &filled_vectors,
                             &more_vectors);

    print_process_duration(better_dot_time, i, la_vector_name,
                           "compute_better_dot_products");


    /*auto [even_better_dot_products,
          even_better_dot_time] = time_this_process<double>
                                  (compute_even_better_dot_products,
                                  &filled_vectors,
                                  &more_vectors);
    std::cout << "Dot via std::transform_reduce time of "
              << i
              << " 4d vectors took "
              << even_better_dot_time
              << std::endl;
    std::cout << "\ti.e., "
              << even_better_dot_time.count() / filled_vectors.size()
              << "ns per vector pair"
              << std::endl;*/
    std::cout << std::endl;
  }
  std::cout << "**\t^ABOVE^\t**\tla_vector" << std::endl;
  std::cout << "**\tBELOW\t**\tstd_vector" << std::endl << std::endl;
  for (int i = 1; i < 100000000; i = i * 10)
  {
    auto [vectors,
          construction_time] = time_this_process<std_vector<double>>
                               (construct_n_std_vectors_of_dim, i, 4);
    print_process_duration(construction_time, i, std_vector_name,
                           "construct_n_std_vectors_of_dim");
    auto [more_vectors,
          list_construction_time] = time_this_process<std_vector<double>>
                                    (list_construct_n_std_vectors, i);
    print_process_duration(list_construction_time, i, std_vector_name,
                           "list_construct_n_vectors");

    auto [filled_vectors,
          fill_time] = time_this_process<std_vector<double>>
                       (fill_zero_std_vectors, std::move(vectors));
    print_process_duration(fill_time, i, std_vector_name,
                           "fill_zero_std_vectors");

    auto [dot_products,
          dot_time] = time_this_process<double>(compute_std_dot_products,
                                        &filled_vectors, &more_vectors);
    print_process_duration(dot_time, i, std_vector_name,
                           "compute_std_dot_products");

    auto [better_dot_products,
          better_dot_time] = time_this_process<double>
                             (compute_better_std_dot_products,
                             &filled_vectors,
                             &more_vectors);
    print_process_duration(better_dot_time, i, std_vector_name,
                          "compute_better_std_dot_products \
                          (std::transform_reduce)");

    std::cout << std::endl;
  }
  std::cout << "END" << std::endl;
}
