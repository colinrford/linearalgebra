option(LAM_ENABLE_ASAN "Enable AddressSanitizer" OFF)
option(LAM_ENABLE_UBSAN "Enable UndefinedBehaviorSanitizer" OFF)
option(LAM_ENABLE_TSAN "Enable ThreadSanitizer (cannot combine with ASan)" OFF)

if(LAM_ENABLE_TSAN AND LAM_ENABLE_ASAN)
  message(FATAL_ERROR "Cannot enable both TSan and ASan simultaneously")
endif()

if(LAM_ENABLE_ASAN)
  add_compile_options(-fsanitize=address -fno-omit-frame-pointer)
  add_link_options(-fsanitize=address)
endif()

if(LAM_ENABLE_UBSAN)
  add_compile_options(-fsanitize=undefined -fno-omit-frame-pointer)
  add_link_options(-fsanitize=undefined)
endif()

if(LAM_ENABLE_TSAN)
  add_compile_options(-fsanitize=thread -fno-omit-frame-pointer)
  add_link_options(-fsanitize=thread)
endif()
