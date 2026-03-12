set(_lib "fortbridge")
set(_pkg "FORTBRIDGE")
set(_url "https://github.com/pprcht/fortbridge")

# Discovery method order can be overridden by the parent project, e.g.:
#   set(FORTBRIDGE_FIND_METHOD "subproject" "cmake")
if(NOT DEFINED "${_pkg}_FIND_METHOD")
  set("${_pkg}_FIND_METHOD" "cmake" "subproject" "fetch" "pkgconf")
endif()

# Reuse whichever utils macro your main project already provides.
# Replace "crest-utils" with the actual name if yours differs.
include("${CMAKE_CURRENT_LIST_DIR}/crest-utils.cmake")

crest_find_package("${_lib}" "${${_pkg}_FIND_METHOD}" "${_url}")

set(found FALSE)
if(TARGET "fortbridge::fortbridge")
  set(found TRUE)
endif()
message(STATUS "Found fortbridge: ${found}")

unset(_lib)
unset(_pkg)
unset(_url)
