AC_INIT([kwip], [1.0], [kdmfoss@gmail.com])
AM_INIT_AUTOMAKE([-Wall foreign])
AC_PROG_CC
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
 Makefile
 src/Makefile
])

PKG_CHECK_MODULES([hdf5], [hdf5])
PKG_CHECK_MODULES([zlib], [zlib])

AC_OUTPUT
