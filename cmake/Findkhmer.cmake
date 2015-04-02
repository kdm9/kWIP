# - Find khmer
# Find the native khmer includes and library.
# Once done this will define
#
#  KHMER_INCLUDE_DIRS   - where to find qes.h, etc.
#  KHMER_LIBRARIES      - List of libraries when using khmer.
#  KHMER_FOUND          - True if khmer found.
#
#  KHMER_VERSION_STRING - The version of khmer found (x.y.z)
#  KHMER_VERSION_MAJOR  - The major version of khmer
#  KHMER_VERSION_MINOR  - The minor version of khmer
#  KHMER_VERSION_PATCH  - The patch version of khmer
#  KHMER_VERSION_PREREL - The pre-release version of khmer
#  KHMER_VERSION_GIT    - The git version of khmer
#
# An includer may set KHMER_ROOT to a khmer installation root to tell
# this module where to look.

set(_KHMER_SEARCHES)

# Search KHMER_ROOT first if it is set.
if(KHMER_ROOT)
  set(_KHMER_SEARCH_ROOT PATHS ${KHMER_ROOT} NO_DEFAULT_PATH)
  list(APPEND _KHMER_SEARCHES _KHMER_SEARCH_ROOT)
endif()

# Normal search.
set(_KHMER_SEARCH_NORMAL
  PATHS "$ENV{PROGRAMFILES}/khmer"
  )
list(APPEND _KHMER_SEARCHES _KHMER_SEARCH_NORMAL)

# Try each search configuration.
foreach(search ${_KHMER_SEARCHES})
  find_path(KHMER_INCLUDE_DIR NAMES khmer.hh ${${search}} PATH_SUFFIXES include)
  find_library(KHMER_LIBRARY NAMES libkhmer.a ${${search}} PATH_SUFFIXES lib)
endforeach()

mark_as_advanced(KHMER_LIBRARY KHMER_INCLUDE_DIR)

set(KHMER_VERSION_MAJOR "1")
set(KHMER_VERSION_MINOR "3")
set(KHMER_VERSION_STRING "${KHMER_VERSION_MAJOR}.${KHMER_VERSION_MINOR}")

# handle the QUIETLY and REQUIRED arguments and set KHMER_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(KHMER REQUIRED_VARS KHMER_LIBRARY KHMER_INCLUDE_DIR
                                       VERSION_VAR KHMER_VERSION_STRING)

if(KHMER_FOUND)
    set(KHMER_INCLUDE_DIRS ${KHMER_INCLUDE_DIR})
    set(KHMER_LIBRARIES ${KHMER_LIBRARY})
endif()


