AC_INIT([kwip], [1.0], [kdmfoss@gmail.com])
AM_INIT_AUTOMAKE([-Wall foreign subdir-objects])
AC_CONFIG_MACRO_DIR([.autoconf_archive/m4])

AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX

PKG_CHECK_MODULES([hdf5], [hdf5])
PKG_CHECK_MODULES([gsl], [gsl])
PKG_CHECK_MODULES([zlib], [zlib])

AX_PTHREAD
AX_EXT

AC_ARG_WITH([sse],
            [AS_HELP_STRING([--with-sse], [enable SSE])],
            [], [with_sse=yes])
AM_CONDITIONAL([HAVE_SSE2], test  "x$with_sse" != xno)
AC_ARG_WITH([avx],
            [AS_HELP_STRING([--with-avx], [enable avx])],
            [], [with_avx=no])
AM_CONDITIONAL([HAVE_AVX2], test  "x$with_avx" != xno)


AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_FILES([
 Makefile
 src/Makefile
])

AC_OUTPUT
