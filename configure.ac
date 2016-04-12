dnl Process this file with autoconf to produce a configure script.
dnl Author: Erik.Bartos@gmail.com

AC_PREREQ([2.69])
AC_INIT([buenavista], [0.4], [Erik.Bartos@gmail.com])
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])

AM_INIT_AUTOMAKE([subdir-objects])
LT_INIT
AM_SILENT_RULES([yes])

# Checks for programs.
AC_PROG_CXX
AC_PROG_AWK
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_PROG_MKDIR_P

# Checks for libraries.

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([sqrt floor])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h unistd.h])

GETTEXT_PACKAGE=buenavista
AC_SUBST(GETTEXT_PACKAGE)
AC_DEFINE_UNQUOTED(GETTEXT_PACKAGE,"$GETTEXT_PACKAGE", [GETTEXT package name])
AM_GLIB_GNU_GETTEXT

AC_CONFIG_FILES([Data/Makefile
                 Docs/Makefile
                 Hyperons/Makefile
                 Makefile
                 Nucleons/Delta/Makefile
                 Nucleons/Iachello/Makefile
                 Nucleons/Makefile
                 Pions/Makefile
                 Transition/Makefile
                 src/Makefile])
AC_OUTPUT

echo
echo " BuenaVista $VERSION"
echo

