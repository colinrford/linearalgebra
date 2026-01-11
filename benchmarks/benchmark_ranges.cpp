import lam.linearalgebra;
import std;

using namespace lam::linalg;
using namespace std::chrono;

// Prevent optimization
volatile double sink;

void bench_vector_construct(int iterations)
{
  constexpr std::size_t vec_size = 10000;
  std::vector<double> data(vec_size);
  std::generate(data.begin(), data.end(), std::rand);

  auto start = steady_clock::now();
  for (int i = 0; i < iterations; ++i)
  {
    vector<double> v(data);
    sink = v[0];
  }
  auto end = steady_clock::now();
  auto dur = duration_cast<microseconds>(end - start).count();
  std::println("vector construct:   {} us", dur);
}

void bench_r_vector_construct(int iterations)
{
  constexpr std::size_t vec_size = 10000;
  std::vector<double> data(vec_size);
  std::generate(data.begin(), data.end(), std::rand);

  auto start = steady_clock::now();
  for (int i = 0; i < iterations; ++i)
  {
    r_vector<double> v(data);
    sink = v[0];
  }
  auto end = steady_clock::now();
  auto dur = duration_cast<microseconds>(end - start).count();
  std::println("r_vector construct: {} us", dur);
}

void bench_vector_add(int iterations)
{
  constexpr std::size_t vec_size = 10000;
  std::vector<double> data(vec_size, 1.0);
  vector<double> v1(data);
  vector<double> v2(data);

  auto start = steady_clock::now();
  for (int i = 0; i < iterations; ++i)
  {
    v1 += v2; // In place add
    sink = v1[0];
  }
  auto end = steady_clock::now();
  auto dur = duration_cast<microseconds>(end - start).count();
  std::println("vector add:         {} us", dur);
}

void bench_r_vector_add(int iterations)
{
  constexpr std::size_t vec_size = 10000;
  std::vector<double> data(vec_size, 1.0);
  r_vector<double> v1(data);
  r_vector<double> v2(data);

  auto start = steady_clock::now();
  for (int i = 0; i < iterations; ++i)
  {
    v1 += v2; // In place add
    sink = v1[0];
  }
  auto end = steady_clock::now();
  auto dur = duration_cast<microseconds>(end - start).count();
  std::println("r_vector add:       {} us", dur);
}

// PMR Benchmarks
void bench_vector_pmr_construct(int iterations)
{
  constexpr std::size_t vec_size = 10000;
  std::vector<double> data(vec_size);
  std::generate(data.begin(), data.end(), std::rand);

  // Create a buffer large enough for all iterations to avoid upstream alloc noise
  // 10k doubles * 8 bytes * 1000 iterations = ~80MB
  std::vector<std::byte> buffer(100 * 1024 * 1024);
  std::pmr::monotonic_buffer_resource pool{buffer.data(), buffer.size()};
  std::pmr::polymorphic_allocator<double> alloc{&pool};

  using pmr_vec = vector<double, std::pmr::polymorphic_allocator<double>>;

  auto start = steady_clock::now();
  for (int i = 0; i < iterations; ++i)
  {
    // We reuse the pool without resetting for simplicity (monotonic)
    // This measures raw "pointer bump" allocation speed + construction cost
    pmr_vec v(data, alloc);
    sink = v[0];
  }
  auto end = steady_clock::now();
  auto dur = duration_cast<microseconds>(end - start).count();
  std::println("vector (PMR) construct:   {} us", dur);
}

void bench_r_vector_pmr_construct(int iterations)
{
  constexpr std::size_t vec_size = 10000;
  std::vector<double> data(vec_size);
  std::generate(data.begin(), data.end(), std::rand);

  std::vector<std::byte> buffer(100 * 1024 * 1024);
  std::pmr::monotonic_buffer_resource pool{buffer.data(), buffer.size()};
  std::pmr::polymorphic_allocator<double> alloc{&pool};

  using pmr_r_vec = r_vector<double, std::pmr::polymorphic_allocator<double>>;

  auto start = steady_clock::now();
  for (int i = 0; i < iterations; ++i)
  {
    pmr_r_vec v(data, alloc);
    sink = v[0];
  }
  auto end = steady_clock::now();
  auto dur = duration_cast<microseconds>(end - start).count();
  std::println("r_vector (PMR) construct: {} us", dur);
}

int main()
{
  std::srand(123);
  int iterations = 1000;

  std::println("Running Ranges Benchmark ({} iterations, 10k elements)", iterations);

  // Default Allocator
  bench_vector_construct(iterations);
  bench_r_vector_construct(iterations);

  // PMR Allocator
  bench_vector_pmr_construct(iterations);
  bench_r_vector_pmr_construct(iterations);

  bench_vector_add(iterations * 10);
  bench_r_vector_add(iterations * 10);

  return 0;
}
