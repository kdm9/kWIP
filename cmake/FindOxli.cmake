# - Find oxli
# Find the native oxli includes and library.
# Once done this will define
#
#  OXLI_INCLUDE_DIRS   - where to find qes.h, etc.
#  OXLI_LIBRARIES      - List of libraries when using oxli.
#  OXLI_FOUND          - True if oxli found.
#
#  OXLI_VERSION_STRING - The version of oxli found (x.y.z)
#  OXLI_VERSION_MAJOR  - The major version of oxli
#  OXLI_VERSION_MINOR  - The minor version of oxli
#  OXLI_VERSION_PATCH  - The patch version of oxli
#  OXLI_VERSION_PREREL - The pre-release version of oxli
#  OXLI_VERSION_GIT    - The git version of oxli
#
# An includer may set OXLI_ROOT to a oxli installation root to tell
# this module where to look.

set(_OXLI_SEARCHES)

# Search OXLI_ROOT first if it is set.
if(OXLI_ROOT)
  set(_OXLI_SEARCH_ROOT PATHS ${OXLI_ROOT} NO_DEFAULT_PATH)
  list(APPEND _OXLI_SEARCHES _OXLI_SEARCH_ROOT)
endif()

# Normal search.
set(_OXLI_SEARCH_NORMAL
  PATHS "$ENV{PROGRAMFILES}/oxli"
        "$ENV{HOME}"
	"/usr/local"
	"/usr"
  )
list(APPEND _OXLI_SEARCHES _OXLI_SEARCH_NORMAL)

# Try each search configuration.
foreach(search ${_OXLI_SEARCHES})
  find_path(OXLI_INCLUDE_DIR NAMES oxli/counting.hh ${${search}} PATH_SUFFIXES include)
  find_library(OXLI_LIBRARY NAMES liboxli.a ${${search}} PATH_SUFFIXES lib)
endforeach()

mark_as_advanced(OXLI_LIBRARY OXLI_INCLUDE_DIR)

set(OXLI_VERSION_MAJOR "2")
set(OXLI_VERSION_MINOR "0")
set(OXLI_VERSION_PATCH "0")
set(OXLI_VERSION_STRING "${OXLI_VERSION_MAJOR}.${OXLI_VERSION_MINOR}.${OXLI_VERSION_PATCH}")

# handle the QUIETLY and REQUIRED arguments and set OXLI_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(
	OXLI
	FOUND_VAR OXLI_FOUND
	REQUIRED_VARS OXLI_LIBRARY OXLI_INCLUDE_DIR
        VERSION_VAR OXLI_VERSION_STRING)

if(OXLI_FOUND)
    set(OXLI_INCLUDE_DIRS ${OXLI_INCLUDE_DIR})
    set(OXLI_LIBRARIES ${OXLI_LIBRARY})
endif()
