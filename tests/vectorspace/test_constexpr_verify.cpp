import std;
import lam.linearalgebra;

int main() {
    // Test constexpr construction and indexing
    constexpr auto test_constexpr = []() {
        lam::linalg::vector<double> v{1.0, 2.0, 3.0};
        return v[1];
    };
    static_assert(test_constexpr() == 2.0);

    // Test constexpr dot product
    constexpr auto test_dot = []() {
        lam::linalg::vector<double> v1{1.0, 2.0};
        lam::linalg::vector<double> v2{3.0, 4.0};
        return v1.dot(v2);
    };
    static_assert(test_dot() == 11.0);

    // Test constexpr norm2
    constexpr auto test_norm2 = []() {
        lam::linalg::vector<double> v{3.0, 4.0};
        return v.norm2();
    };
    static_assert(test_norm2() == 25.0);

    // Test constexpr norm (sqrt)
    constexpr auto test_norm = []() {
        lam::linalg::vector<double> v{3.0, 4.0};
        return v.norm();
    };
    static_assert(test_norm() == 5.0);

    std::println("Constexpr tests passed!");
    return 0;
}
