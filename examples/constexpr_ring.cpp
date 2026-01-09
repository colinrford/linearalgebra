import std;
import lam.linearalgebra;

using namespace lam::linalg;

// A simple fixed-size ring buffer resource
template<std::size_t Size>
class RingResource
{
  std::array<std::byte, Size> buffer;
  std::size_t offset = 0;

public:
  constexpr RingResource() = default;

  // Not constexpr because reinterpret_cast is invalid in constexpr
  void* allocate(std::size_t bytes, std::size_t align)
  {
    // Simple ring logic: bump offset, wrap if full?
    // For simplicity: Linear Monotonic that resets.
    std::size_t padding = 0;
    void* ptr = buffer.data() + offset;
    std::size_t space = Size - offset;
    // Use std::align to adjust ptr?
    if (std::align(align, bytes, ptr, space))
    {
      offset = ((std::byte*)ptr - buffer.data()) + bytes;
      return ptr;
    }
    // Out of memory: Reset to beginning (Ring behavior)
    offset = 0;
    ptr = buffer.data();
    space = Size;
    if (std::align(align, bytes, ptr, space))
    {
      offset = ((std::byte*)ptr - buffer.data()) + bytes;
      return ptr;
    }
    throw std::bad_alloc();
  }
  void deallocate(void*, std::size_t, std::size_t)
  {
    // No-op
  }
};

template<typename T, std::size_t Size>
struct RingAllocator
{
  using value_type = T;

  RingResource<Size>* resource = nullptr;

  constexpr RingAllocator() = default;
  constexpr RingAllocator(RingResource<Size>& r) : resource(&r) {}

  template<typename U>
  constexpr RingAllocator(const RingAllocator<U, Size>& other) noexcept : resource(other.resource)
  {}

  constexpr T* allocate(std::size_t n)
  {
    if (std::is_constant_evaluated())
    {
      return std::allocator<T>{}.allocate(n);
    }
    else
    {
      if (!resource)
        throw std::bad_alloc();
      return static_cast<T*>(resource->allocate(n * sizeof(T), alignof(T)));
    }
  }

  constexpr void deallocate(T* p, std::size_t n)
  {
    if (std::is_constant_evaluated())
    {
      std::allocator<T>{}.deallocate(p, n);
    }
    else
    {
      resource->deallocate(p, n * sizeof(T), alignof(T));
    }
  }

  // Equality
  friend constexpr bool operator==(const RingAllocator& a, const RingAllocator& b)
  {
    if (std::is_constant_evaluated())
      return true; // Stateless at compile time
    return a.resource == b.resource;
  }
};

int main()
{
  // 1. Runtime Usage with Ring Buffer
  {
    RingResource<1024> ring;
    RingAllocator<double, 1024> alloc(ring);

    vector<double, RingAllocator<double, 1024>> v1({1.0, 2.0, 3.0}, alloc);
    vector<double, RingAllocator<double, 1024>> v2({4.0, 5.0, 6.0}, alloc);

    auto v3 = v1 + v2; // Propagates alloc

    std::println("Runtime Ring Calculation: {}", v3[0]);

    if (v3.get_allocator().resource != &ring)
    {
      std::println("Error: Allocator not propagated!");
      return 1;
    }
  }

  // 2. Constexpr Usage (Falls back to heap/std::allocator)
  static_assert([]() constexpr {
    // Can't use RingResource at compile time easily due to reinterpret_cast/buffer limitation
    // So we default construct allocator (resource=nullptr).
    // allocate() checks is_constant_evaluated -> uses std::allocator.
    RingAllocator<double, 1024> alloc; // Null resource

    vector<double, RingAllocator<double, 1024>> v1({10.0, 20.0}, alloc);
    vector<double, RingAllocator<double, 1024>> v2({30.0, 40.0}, alloc);

    auto v3 = v1 + v2;
    return v3[0];
  }() == 40.0);

  std::println("Constexpr Validation Passed!");
  return 0;
}
