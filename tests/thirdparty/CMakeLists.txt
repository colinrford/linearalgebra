# Tests that need a third-party allocator. Nothing here is required: the whole
# directory is skipped when the package is absent, so `tests/` proper stays
# dependency-free (see the minimal-dependency ethos; same pattern as
# benchmarks/micro with google/benchmark).
#
# These are correctness tests -- does lam.linalg behave when handed someone
# else's allocator -- not performance comparisons. Allocator benchmarks belong
# in benchmarks/.

find_package(TBB QUIET)
find_package(Boost QUIET CONFIG COMPONENTS container)


# Shared scaffolding MODULE for the third-party tests: std + lam.linearalgebra only, no Boost
# and no TBB. Built only when something below will use it.
if(TBB_FOUND OR Boost_FOUND)
  add_library(linalg_test_pmr_suite STATIC)
  target_sources(linalg_test_pmr_suite
    PUBLIC FILE_SET CXX_MODULES FILES pmr_resource_suite.cppm)
  target_link_libraries(linalg_test_pmr_suite
    PUBLIC lam_linearalgebra::linearalgebra lam_concepts::concepts)
  target_compile_features(linalg_test_pmr_suite PUBLIC cxx_std_23)
endif()

if(TBB_FOUND)
  # Named for what it is: TBB here supplies allocators to test against, and is NOT linked into
  # the library, nor is it the TBB parallelism backend (that is LAM_USE_TBB, default OFF).
  message(STATUS "linalg: TBB ${TBB_VERSION} found -- building TBB allocator tests (tests only, lib not linked)")
  linalg_add_test(linalg_tbb_allocators_test
    SOURCE test_tbb_allocators.cpp
    LIBS   TBB::tbb TBB::tbbmalloc linalg_test_pmr_suite)
else()
  message(STATUS "linalg: TBB not found -- skipping TBB allocator tests")
endif()

if(Boost_FOUND)
  message(STATUS "linalg: Boost ${Boost_VERSION} found -- building Boost pmr/allocator tests (tests only, lib not linked)")
  linalg_add_test(linalg_boost_pmr_test
    SOURCE test_boost_pmr.cpp
    LIBS   Boost::container linalg_test_pmr_suite)
else()
  message(STATUS "linalg: Boost not found -- skipping Boost pmr/allocator tests")
endif()
