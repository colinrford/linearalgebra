# LamConstexprLimit.cmake
# Portable constexpr evaluation limit: Clang uses -fconstexpr-steps,
# GCC uses -fconstexpr-ops-limit.

function(lam_target_constexpr_limit target limit)
  target_compile_options(${target} PRIVATE
    "$<$<CXX_COMPILER_ID:Clang,AppleClang>:-fconstexpr-steps=${limit}>"
    "$<$<CXX_COMPILER_ID:GNU>:-fconstexpr-ops-limit=${limit}>"
  )
endfunction()
