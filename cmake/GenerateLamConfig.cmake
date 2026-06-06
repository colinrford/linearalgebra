# Build-time generator for <name>_config.cppm.
#
# Invoked by add_custom_command in the submodule's CMakeLists.txt. Re-resolves
# git_hash from HEAD on every invocation, then runs configure_file() against
# the template. All other @VARS@ are passed in by the add_custom_command via
# -D on the cmake -P command line.

set(GIT_HASH "unknown")
if(EXISTS "${SOURCE_DIR}/.git")
  execute_process(
    COMMAND git -C "${SOURCE_DIR}" rev-parse --short HEAD
    OUTPUT_VARIABLE GIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
endif()

configure_file("${TEMPLATE}" "${OUTPUT}" @ONLY)
