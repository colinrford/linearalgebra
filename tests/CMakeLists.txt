
# Example/Test executable for initial verification
add_executable(la_test 
  vectorspace/test_pnorm.cpp
)
target_link_libraries(la_test PRIVATE linearalgebra lam::concepts)

add_executable(la_constexpr_test 
  vectorspace/test_constexpr_verify.cpp
)
target_link_libraries(la_constexpr_test PRIVATE linearalgebra lam::concepts)

add_executable(la_concepts_test 
  transformations/test_matrix_concepts.cpp
)
target_link_libraries(la_concepts_test PRIVATE linearalgebra lam::concepts)

add_executable(la_vectorspace_concepts_test
  vectorspace/test_concepts_verify.cpp
)
target_link_libraries(la_vectorspace_concepts_test PRIVATE linearalgebra lam::concepts)

add_executable(la_test_constexpr_blas
  transformations/test_constexpr_blas.cpp
)
target_link_libraries(la_test_constexpr_blas PRIVATE linearalgebra lam::concepts)

add_executable(test_foreign_matrix 
  transformations/test_foreign_matrix.cpp
)
target_link_libraries(test_foreign_matrix PRIVATE linearalgebra lam::concepts)

add_executable(test_matrix_layout
  transformations/test_matrix_layout.cpp
)
target_link_libraries(test_matrix_layout PRIVATE linearalgebra lam::concepts)


add_executable(la_math_ops_test vectorspace/test_math_ops.cpp)
target_link_libraries(la_math_ops_test PRIVATE linearalgebra lam::concepts)

add_executable(la_allocator_test 
  vectorspace/test_allocator_propagation.cpp
)
target_link_libraries(la_allocator_test PRIVATE linearalgebra lam::concepts)
