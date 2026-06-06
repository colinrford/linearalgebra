# Unit/verification tests. linalg_add_test (cmake/LinalgHelpers.cmake) handles the
# lam::linearalgebra + lam::concepts link and C++23, then registers a ctest.

linalg_add_test(linalg_test                       SOURCE vectorspace/test_pnorm.cpp)
linalg_add_test(linalg_constexpr_test             SOURCE vectorspace/test_constexpr_verify.cpp)
linalg_add_test(linalg_concepts_test              SOURCE transformations/test_matrix_concepts.cpp)
linalg_add_test(linalg_vectorspace_concepts_test  SOURCE vectorspace/test_concepts_verify.cpp)
linalg_add_test(linalg_test_constexpr_blas        SOURCE transformations/test_constexpr_blas.cpp)
linalg_add_test(test_foreign_matrix               SOURCE transformations/test_foreign_matrix.cpp)
linalg_add_test(test_matrix_layout                SOURCE transformations/test_matrix_layout.cpp)
linalg_add_test(linalg_math_ops_test              SOURCE vectorspace/test_math_ops.cpp)
linalg_add_test(linalg_allocator_test             SOURCE vectorspace/test_allocator_propagation.cpp)
linalg_add_test(linalg_test_blas_dispatch         SOURCE vectorspace/test_blas_dispatch.cpp)
linalg_add_test(test_matrix_svd                   SOURCE transformations/test_matrix_svd.cpp)
linalg_add_test(test_matrix_svd_stress            SOURCE transformations/test_matrix_svd_stress.cpp)
linalg_add_test(linalg_test_exceptions            SOURCE vectorspace/test_exceptions.cpp)
linalg_add_test(test_matrix_product_fallback      SOURCE transformations/test_matrix_product_fallback.cpp)
