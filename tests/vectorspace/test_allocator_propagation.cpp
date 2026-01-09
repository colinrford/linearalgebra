import std;
import lam.linearalgebra;
import lam.concepts;

// Helper to check conditions
void check(bool condition, std::string_view msg) {
    if (!condition) {
        throw std::runtime_error(std::string(msg));
    }
}

int main() {
    using namespace lam::linalg;
    
    // Create a buffer for PMR
    std::array<std::byte, 4096> buffer;
    std::pmr::monotonic_buffer_resource pool{buffer.data(), buffer.size()};
    std::pmr::polymorphic_allocator<double> alloc{&pool};
    
    using vec_t = vector<double, std::pmr::polymorphic_allocator<double>>;
    
    // Test 1: Construction with allocator
    {
        vec_t v({1.0, 2.0, 3.0}, alloc);
        check(v.get_allocator() == alloc, "Construction: Allocator not stored correctly");
        std::println("Construction passed");
    }
    
    // Test 2: Copy propagation
    {
        vec_t v1({1.0, 2.0, 3.0}, alloc);
        vec_t v2 = v1; // Copy constructor
        check(v2.get_allocator() == alloc, "Copy: Allocator not propagated");
        std::println("Copy propagation passed");
    }
    
    // Test 3: Allocator-extended copy
    {
        vec_t v1({1.0, 2.0, 3.0}, alloc);
        // Create a different allocator (null resource for testing difference)
        std::pmr::polymorphic_allocator<double> alloc2{}; 
        vec_t v2(v1, alloc2);
        check(v2.get_allocator() == alloc2, "Extended copy: Allocator not used");
        check(v2.get_allocator() != v1.get_allocator(), "Extended copy: Allocator should be different");
        std::println("Allocator-extended copy passed");
    }
    
    // Test 4: Operator+ propagation
    {
        vec_t v1({1.0, 2.0, 3.0}, alloc);
        vec_t v2({4.0, 5.0, 6.0}, alloc);
        vec_t v3 = v1 + v2;
        check(v3.get_allocator() == alloc, "Operator+: Allocator not propagated");
        std::println("Operator+ passed");
    }
    
    // Test 5: Operator* propagation
    {
        vec_t v1({1.0, 2.0, 3.0}, alloc);
        vec_t v2 = v1 * 2.0;
        check(v2.get_allocator() == alloc, "Operator*: Allocator not propagated");
        vec_t v3 = 2.0 * v1;
        check(v3.get_allocator() == alloc, "Operator* (scalar first): Allocator not propagated");
        std::println("Operator* passed");
    }
    
    // Test 6: Member functions (cross, unit, etc.)
    {
        vec_t v1({1.0, 0.0, 0.0}, alloc);
        vec_t v2({0.0, 1.0, 0.0}, alloc);
        
        vec_t v_cross = v1.cross(v2);
        check(v_cross.get_allocator() == alloc, "cross(): Allocator not propagated");
        
        vec_t v_unit = v1.unit();
        check(v_unit.get_allocator() == alloc, "unit(): Allocator not propagated");
        
        vec_t v_project = v1.project(v2);
        check(v_project.get_allocator() == alloc, "project(): Allocator not propagated"); // Should take from *this (v1)
        
        vec_t v_lerp = v1.lerp(v2, 0.5);
        check(v_lerp.get_allocator() == alloc, "lerp(): Allocator not propagated");
        
        std::println("Member functions passed");
    }

    // Test 7: Propagation vs Non-Propagation (Comparison)
    {
        vec_t v1({10.0, 20.0}, alloc);

        // Case A: Propagation (Default Copy Constructor)
        vec_t v_prop = v1;
        check(v_prop.get_allocator() == alloc, "Case A: Default copy should propagate allocator");
        check(v_prop.get_allocator() == v1.get_allocator(), "Case A: Allocators should match");

        // Case B: Non-Propagation (Explicitly using default allocator)
        std::pmr::polymorphic_allocator<double> default_alloc{}; // Uses default resource (new/delete)
        vec_t v_noprop(v1, default_alloc);
        
        check(v_noprop.get_allocator() == default_alloc, "Case B: Explicit copy should use provided allocator");
        check(v_noprop.get_allocator() != alloc, "Case B: Should NOT use source allocator");
        check(v_noprop.get_allocator().resource() != alloc.resource(), "Case B: Resources should differ (Default vs Monotonic)");

        std::println("Propagation vs Non-Propagation verified");
    }

    std::println("\n✅ All allocator propagation tests passed!");
    return 0;
}
