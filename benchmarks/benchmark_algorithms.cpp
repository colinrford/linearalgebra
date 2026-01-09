
import std;
import lam.linearalgebra;

using namespace lam::linalg;

template<typename Func>
auto measure(std::string_view name, Func&& func)
{
    auto start = std::chrono::high_resolution_clock::now();
    auto result = func();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> ms = end - start;
    std::println("{:<30} {:10.4f} ms (result: {})", name, ms.count(), result);
    return ms.count();
}

int main()
{
    constexpr std::size_t N = 10'000'000;
    std::println("Benchmarking with N = {}", N);

    // Setup data using std::vector for raw data source
    std::vector<double> v1_data(N, 1.0);
    std::vector<double> v2_data(N, 2.0);

    // Create lam::linalg::vector (legacy) - copies data
    vector<double> vec1(v1_data);
    vector<double> vec2(v2_data);

    // Create lam::linalg::r_vector (range) - copies data
    r_vector<double> rvec1(v1_data);
    r_vector<double> rvec2(v2_data);

    std::println("\n--- Dot Product ---");
    measure("vector::dot (legacy)", [&]() { return vec1.dot(vec2); });
    measure("r_vector::dot (legacy)", [&]() { return rvec1.dot(rvec2); });
    measure("dot_range(vector)", [&]() { return dot_range(vec1, vec2); });
    measure("dot_range(r_vector)", [&]() { return dot_range(rvec1, rvec2); });
    measure("dot_range(std::vector)", [&]() { return dot_range(v1_data, v2_data); });

    std::println("\n--- Distance ---");
    measure("vector::distance (legacy)", [&]() { return vec1.distance(vec2); });
    measure("r_vector::distance (legacy)", [&]() { return rvec1.distance(rvec2); });
    measure("distance_range(vector)", [&]() { return distance_range(vec1, vec2); });
    measure("distance_range(r_vector)", [&]() { return distance_range(rvec1, rvec2); });

    std::println("\n--- Norm ---");
    measure("vector::norm (legacy)", [&]() { return vec1.norm(); });
    measure("r_vector::norm (legacy)", [&]() { return rvec1.norm(); });
    measure("norm_range(vector)", [&]() { return norm_range(vec1); });
    measure("norm_range(r_vector)", [&]() { return norm_range(rvec1); });
    
    return 0;
}
