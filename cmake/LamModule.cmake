# LamModule.cmake — shared lam submodule build helpers.
#
# Master copy lives at <workspace>/cmake/LamModule.cmake. Each submodule carries
# its OWN vendored copy in <submodule>/cmake/ so a submodule still configures
# standalone with no cross-submodule path dependency (same convention as
# GenerateLamConfig.cmake / Sanitizers.cmake). Keep the copies in sync from the
# master; do not include the master via an out-of-tree path.
#
# Provides:
#   lam_add_config_partition(<target> NAME <name> [TEMPLATE <file>]
#                            [EXTRA_DEFINES <-Dk=v>...])
#   lam_install_module_package(<target> NAME <name>
#                              DESCRIPTION <text> LICENSE <spdx>
#                              [MODULE_SUBDIR <dir>] [REQUIRES <pkg>...])
#   lam_add_config_dump(<link-target> NAME <name>)
#   lam_find_dependency(<pkg> [COMPONENT <comp>])

if(_LAM_MODULE_INCLUDED)
  return()
endif()
set(_LAM_MODULE_INCLUDED TRUE)

# -----------------------------------------------------------------------------
# lam_add_config_partition(<target> NAME <name> [TEMPLATE <abs-or-rel-path>]
#                          [EXTRA_DEFINES <-Dkey=val>...])
#
# Generates the lam.<name>:config module partition (<name>_config.cppm) at BUILD
# time and appends it to <target>'s CXX_MODULES file set.
#
# The generated file's git_hash refreshes whenever .git/logs/HEAD changes
# (commit / checkout / reset / amend ...) WITHOUT requiring a reconfigure: the
# add_custom_command depends on logs/HEAD and re-runs GenerateLamConfig.cmake at
# build time. All other fields are configure-time substitutions.
#
# EXTRA_DEFINES passes additional `-Dkey=val` flags to the worker script for
# submodules whose template has fields beyond the standard set (e.g.
# polynomial_nttp's BLAS / Accelerate / TBB / int128 / TSAN booleans).
#
# Requires a vendored cmake/GenerateLamConfig.cmake worker script. TEMPLATE
# defaults to ${CMAKE_CURRENT_SOURCE_DIR}/<name>_config.cppm.in.
# -----------------------------------------------------------------------------
function(lam_add_config_partition target)
  cmake_parse_arguments(PARSE_ARGV 1 _arg "" "NAME;TEMPLATE" "EXTRA_DEFINES")
  if(NOT _arg_NAME)
    message(FATAL_ERROR "lam_add_config_partition: NAME is required")
  endif()
  if(NOT _arg_TEMPLATE)
    set(_arg_TEMPLATE "${CMAKE_CURRENT_SOURCE_DIR}/${_arg_NAME}_config.cppm.in")
  endif()

  # Resolve the real logs/HEAD path. `git rev-parse --git-path` handles the
  # gitlink case (submodules where .git is a file pointing into ../.git/modules/).
  execute_process(
    COMMAND git -C "${CMAKE_CURRENT_SOURCE_DIR}" rev-parse --git-path logs/HEAD
    OUTPUT_VARIABLE _git_logs_head
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET)
  get_filename_component(_git_logs_head "${_git_logs_head}"
    ABSOLUTE BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

  if(CMAKE_BUILD_TYPE AND CMAKE_BUILD_TYPE MATCHES "^[Dd][Ee][Bb][Uu][Gg]$")
    set(_is_debug "true")
  else()
    set(_is_debug "false")
  endif()

  set(_cfg_out "${CMAKE_CURRENT_BINARY_DIR}/${_arg_NAME}_config.cppm")
  add_custom_command(
    OUTPUT  ${_cfg_out}
    COMMAND ${CMAKE_COMMAND}
            -DTEMPLATE=${_arg_TEMPLATE}
            -DOUTPUT=${_cfg_out}
            -DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
            -DPROJECT_VERSION=${PROJECT_VERSION}
            -DPROJECT_VERSION_MAJOR=${PROJECT_VERSION_MAJOR}
            -DPROJECT_VERSION_MINOR=${PROJECT_VERSION_MINOR}
            -DPROJECT_VERSION_PATCH=${PROJECT_VERSION_PATCH}
            -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
            -DIS_DEBUG=${_is_debug}
            -DCMAKE_CXX_COMPILER_ID=${CMAKE_CXX_COMPILER_ID}
            -DCMAKE_CXX_COMPILER_VERSION=${CMAKE_CXX_COMPILER_VERSION}
            -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
            -DCMAKE_SYSTEM_NAME=${CMAKE_SYSTEM_NAME}
            -DCMAKE_SYSTEM_PROCESSOR=${CMAKE_SYSTEM_PROCESSOR}
            -DCMAKE_GENERATOR=${CMAKE_GENERATOR}
            -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
            "-DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}"
            ${_arg_EXTRA_DEFINES}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/GenerateLamConfig.cmake
    DEPENDS ${_arg_TEMPLATE}
            ${_git_logs_head}
    VERBATIM)

  target_sources(${target} PUBLIC
    FILE_SET CXX_MODULES TYPE CXX_MODULES
    BASE_DIRS ${CMAKE_CURRENT_BINARY_DIR}
    FILES ${_cfg_out})
endfunction()

# -----------------------------------------------------------------------------
# lam_install_module_package(<target> NAME <name>
#                            DESCRIPTION <text> LICENSE <spdx>
#                            [MODULE_SUBDIR <dir>] [REQUIRES <pkg>...])
#
# Installs <target> as package lam_<name> following the Path A convention:
#   - archive/runtime + CXX_MODULES file set -> include/lam/<subdir>
#   - install(EXPORT) with NAMESPACE lam_<name>::  (matches find_package name)
#   - install(PACKAGE_INFO) primary CPS, guarded for CMake >= 4.3
#   - legacy lam_<name>Config.cmake + ConfigVersion.cmake for older tooling
# VERSION / COMPAT_VERSION derive from PROJECT_VERSION (major.minor.0 floor).
# MODULE_SUBDIR defaults to <name>.
# REQUIRES lists transitive lam packages emitted as find_dependency(<pkg>) lines
# in the legacy Config.cmake (e.g. REQUIRES lam_concepts for a consumer module).
# -----------------------------------------------------------------------------
function(lam_install_module_package target)
  cmake_parse_arguments(PARSE_ARGV 1 _arg ""
    "NAME;DESCRIPTION;LICENSE;MODULE_SUBDIR" "REQUIRES")
  if(NOT _arg_NAME)
    message(FATAL_ERROR "lam_install_module_package: NAME is required")
  endif()
  if(NOT _arg_MODULE_SUBDIR)
    set(_arg_MODULE_SUBDIR "${_arg_NAME}")
  endif()

  include(GNUInstallDirs)
  include(CMakePackageConfigHelpers)

  set(_pkg "lam_${_arg_NAME}")
  set(_targets "${_pkg}Targets")
  set(_moddir "${CMAKE_INSTALL_INCLUDEDIR}/lam/${_arg_MODULE_SUBDIR}")

  install(TARGETS ${target}
    EXPORT ${_targets}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    FILE_SET CXX_MODULES DESTINATION ${_moddir})

  install(EXPORT ${_targets}
    FILE ${_targets}.cmake
    NAMESPACE ${_pkg}::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${_pkg})

  # Primary CPS (Common Package Specification). Guarded so CMake < 4.3, which
  # lacks install(PACKAGE_INFO), still configures (legacy *Config.cmake remains).
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 4.3)
    install(PACKAGE_INFO ${_pkg}
      EXPORT       ${_targets}
      VERSION      ${PROJECT_VERSION}
      COMPAT_VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.0
      DESCRIPTION  "${_arg_DESCRIPTION}"
      LICENSE      "${_arg_LICENSE}"
      CXX_MODULES_DIRECTORY "include/lam/${_arg_MODULE_SUBDIR}"
      DESTINATION  ${CMAKE_INSTALL_LIBDIR}/cps)
  endif()

  # Legacy *Config.cmake: pulls transitive lam deps via find_dependency, then
  # includes the exported targets file. \${CMAKE_CURRENT_LIST_DIR} stays literal
  # (escaped) so it resolves at the consumer's find_package time, not now.
  set(_finddeps "")
  foreach(_dep IN LISTS _arg_REQUIRES)
    string(APPEND _finddeps "find_dependency(${_dep})\n")
  endforeach()
  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/${_pkg}Config.cmake"
"include(CMakeFindDependencyMacro)\n${_finddeps}include(\"\${CMAKE_CURRENT_LIST_DIR}/${_targets}.cmake\")\n")

  write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/${_pkg}ConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion)

  install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/${_pkg}Config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/${_pkg}ConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${_pkg})
endfunction()

# -----------------------------------------------------------------------------
# lam_add_config_dump(<link-target> NAME <name>)
#
# Generates a tiny <name>_config_dump executable that imports lam.<name> and
# calls lam::<name>::config::print(), and registers it as a ctest. Doubles as a
# smoke test (the partition must import + link) and a build-config diagnostic.
# The module's :config partition must define a print() (see *_config.cppm.in).
# Requires enable_testing() in the calling submodule. <link-target> is the
# module library target to link (e.g. the bare target or its alias).
# -----------------------------------------------------------------------------
function(lam_add_config_dump link_target)
  cmake_parse_arguments(PARSE_ARGV 1 _arg "" "NAME;MODULE;NAMESPACE" "")
  if(NOT _arg_NAME)
    message(FATAL_ERROR "lam_add_config_dump: NAME is required")
  endif()

  # MODULE/NAMESPACE default to the lam.<name> / lam::<name>::config convention,
  # but may be overridden when the module name and namespace diverge from <name>
  # (e.g. polynomial_nttp: module lam.polynomial.nttp, ns lam::polynomial::nttp).
  if(NOT _arg_MODULE)
    set(_arg_MODULE "lam.${_arg_NAME}")
  endif()
  if(NOT _arg_NAMESPACE)
    set(_arg_NAMESPACE "lam::${_arg_NAME}::config")
  endif()

  set(_src "${CMAKE_CURRENT_BINARY_DIR}/${_arg_NAME}_config_dump.cpp")
  file(WRITE "${_src}"
"import ${_arg_MODULE};\nint main() { ${_arg_NAMESPACE}::print(); return 0; }\n")

  add_executable(${_arg_NAME}_config_dump "${_src}")
  target_link_libraries(${_arg_NAME}_config_dump PRIVATE ${link_target})
  target_compile_features(${_arg_NAME}_config_dump PRIVATE cxx_std_23)
  add_test(NAME ${_arg_NAME}_config_dump COMMAND ${_arg_NAME}_config_dump)
endfunction()

# -----------------------------------------------------------------------------
# lam_find_dependency(<pkg> [COMPONENT <comp>])
#
# find_package(<pkg>) for a lam dependency, then applies the BMI-synthesis
# property dance the consumer needs on the imported target. COMPONENT defaults
# to <pkg> with a leading "lam_" stripped (e.g. lam_concepts -> concepts), so
# the imported target is <pkg>::<comp> (lam_concepts::concepts).
#
# Replaces the per-consumer boilerplate:
#   if(NOT TARGET lam_x::x) find_package(lam_x REQUIRED) endif()
#   set_target_properties(lam_x::x PROPERTIES CXX_MODULE_STD ON CXX_EXTENSIONS OFF)
#
# Without the property dance the synthesized BMI compile of the imported
# module fails with "module 'std' not found" or a -std=gnu++23/c++23 mismatch.
# -----------------------------------------------------------------------------
function(lam_find_dependency pkg)
  cmake_parse_arguments(PARSE_ARGV 1 _arg "" "COMPONENT" "")
  set(_comp "${_arg_COMPONENT}")
  if(NOT _comp)
    string(REGEX REPLACE "^lam_" "" _comp "${pkg}")
  endif()
  set(_tgt "${pkg}::${_comp}")

  if(NOT TARGET ${_tgt})
    find_package(${pkg} REQUIRED)
  endif()
  if(TARGET ${_tgt})
    set_target_properties(${_tgt} PROPERTIES CXX_MODULE_STD ON CXX_EXTENSIONS OFF)
  endif()
endfunction()
