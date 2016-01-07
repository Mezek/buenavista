Idnl Process this file with autoconf to produce a configure script.
dnl Author: Erik.Bartos@gmail.com

AC_INIT(buenavista, 0.3)

AC_CONFIG_HEADERS([config.h])

AM_INIT_AUTOMAKE([1.15])

AM_SILENT_RULES([yes])

AC_PROG_CXX

GETTEXT_PACKAGE=buenavista
AC_SUBST(GETTEXT_PACKAGE)
AC_DEFINE_UNQUOTED(GETTEXT_PACKAGE,"$GETTEXT_PACKAGE", [GETTEXT package name])
AM_GLIB_GNU_GETTEXT

LT_INIT

AC_OUTPUT([
Makefile
Data/Makefile
Docs/Makefile
Hyperons/Makefile
Nucleons/Makefile
Pions/Makefile
src/Makefile
Nucleons/Iachello/Makefile
Nucleons/Delta/Makefile
])

echo
echo " BuenaVista $VERSION"
echo

