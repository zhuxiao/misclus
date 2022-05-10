#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([misEval], [0.1.0], [1668191303@qq.com])
AM_INIT_AUTOMAKE([-Wportability])
AC_CONFIG_SRCDIR([Paras.cpp])
AC_CONFIG_HEADERS([config.h])

AC_CHECK_HEADER(miseval/util.h)
AC_CONFIG_MACRO_DIRS([m4])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC


# Checks for libraries.
# FIXME: Replace `main' with a function in `-lhts':
AC_CHECK_LIB([hts], [main],[],[echo "you need to have libhts around"])
# FIXME: Replace `main' with a function in `-lpthread':
AC_CHECK_LIB([pthread], [main],[],[echo "you need to have libpthread around"])

LT_INIT

AM_INIT_AUTOMAKE
AC_PROG_LIBTOOL
AC_LIBTOOL_DLOPEN
AC_PROG_LIBTOOL


# Checks for header files.
AC_CHECK_HEADERS()
AC_CHECK_FUNCS([pow])
AC_CHECK_FUNCS([sqrt])
AC_TYPE_UINT8_T

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_INT8_T
AC_TYPE_SIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T

# Checks for library functions.
AC_FUNC_LSTAT_FOLLOWS_SLASHED_SYMLINK
AC_FUNC_MALLOC
AC_CHECK_FUNCS([getcwd mkdir])

AC_CONFIG_FILES(Makefile)

AC_OUTPUT
