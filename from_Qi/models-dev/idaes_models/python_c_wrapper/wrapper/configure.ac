dnl Initialize version
AC_INIT(libwrapper, version-0.1)
AC_PREREQ([2.69])

AC_PROG_CXX
AC_PROG_CC_C99

dnl To use with libtool and automake:
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([1.14 -Wall no-define foreign])

dnl On Docker image, use export CPPFLAGS="-I/usr/local/src/Ipopt-3.12.5/ThirdParty/ASL -I/usr/local/src/Ipopt-3.12.5/ThirdParty/ASL/solvers -I/opt/conda/pkgs/python-2.7.12-0/include/python2.7"

dnl Checking for ADOL-C header:
AC_CHECK_HEADERS(["solvers/funcadd.h"],[],
                [AC_MSG_ERROR([
ERROR:  Missing ADOL-C headers, Add the proper include path to CPPFLAGS with '-I'
Fatal:  ADOL-C headers not found.])],
                [])

dnl Check Python header:
AC_CHECK_HEADERS([Python.h], [],
                [AC_MSG_ERROR([
ERROR:  Missing Python headers, Add the proper include path to CPPFLAGS with '-I'
Fatal:  Python headers not found.])])

dnl Set ADOL-C directory linker flag:
AC_ARG_WITH(adolc,
        AC_HELP_STRING(
                [--with-adolc=PROGRAM],
                [location of the ADOL-C library]),
                [AC_SUBST(ADOLC_LIBS,"-L$withval")],
                [AC_SUBST(ADOLC_LIBS,"")])
dnl AC_SUBST(ADOLC_LIBS,"-L$adolc_prog")

dnl Set Python directory linker flag (needed only on NERSC?)
AC_ARG_WITH(python,
        AC_HELP_STRING(
                [--with-python=PROGRAM],
                [location of the Python library]),
                [AC_SUBST(PYTHON_LIBS,"-L$withval")],
                [AC_SUBST(PYTHON_LIBS,"")])
dnl AC_SUBST(PYTHON_LIBS,"-L$python_prog")

dnl Libtool and automake setup:
AM_PROG_AR
LT_INIT
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
