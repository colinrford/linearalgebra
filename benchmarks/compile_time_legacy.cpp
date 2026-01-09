
import std;
import lam.linearalgebra;

// Stress test compile time by instantiating many constexpr calls
constexpr double run_legacy() {
    double sum = 0;
    // Repeat many times to stress constant evaluator
    for (int i = 0; i < 500; ++i) {
        lam::linalg::vector<double> v1{1.0, 2.0, 3.0};
        lam::linalg::vector<double> v2{2.0, 3.0, 4.0};
        sum += lam::linalg::dot(v1, v2);
        sum += lam::linalg::distance(v1, v2);
        sum += lam::linalg::norm(v1);
    }
    return sum;
}

static_assert(run_legacy() > 0);

int main() {
    return 0;
}
