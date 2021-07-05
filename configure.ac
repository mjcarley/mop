AC_INIT([mop], 0.1)
AC_CONFIG_SRCDIR([src/mop.h])

# Making releases:
#   MOP_MICRO_VERSION += 1;
#   MOP_INTERFACE_AGE += 1;
#   MOP_BINARY_AGE += 1;
# if any functions have been added, set MOP_INTERFACE_AGE to 0.
# if backwards compatibility has been broken,
# set MOP_BINARY_AGE and MOP_INTERFACE_AGE to 0.
#
MOP_MAJOR_VERSION=0
MOP_MINOR_VERSION=1
MOP_MICRO_VERSION=0
MOP_INTERFACE_AGE=0
MOP_BINARY_AGE=0
MOP_VERSION=$MOP_MAJOR_VERSION.$MOP_MINOR_VERSION.$MOP_MICRO_VERSION

AC_SUBST(MOP_MAJOR_VERSION)
AC_SUBST(MOP_MINOR_VERSION)
AC_SUBST(MOP_MICRO_VERSION)
AC_SUBST(MOP_VERSION)

# libtool versioning
LT_RELEASE=$MOP_MAJOR_VERSION.$MOP_MINOR_VERSION
LT_CURRENT=`expr $MOP_MICRO_VERSION - $MOP_INTERFACE_AGE`
LT_REVISION=$MOP_INTERFACE_AGE
LT_AGE=`expr $MOP_BINARY_AGE - $MOP_INTERFACE_AGE`
AC_SUBST(LT_RELEASE)
AC_SUBST(LT_CURRENT)
AC_SUBST(LT_REVISION)
AC_SUBST(LT_AGE)

# For automake.
VERSION=$MOP_VERSION
PACKAGE=mop

AM_INIT_AUTOMAKE(subdir-objects)
##$PACKAGE, $VERSION, no-define)

AC_DEFINE_UNQUOTED(MOP_MAJOR_VERSION, $MOP_MAJOR_VERSION, [Major version])
AC_DEFINE_UNQUOTED(MOP_MINOR_VERSION, $MOP_MINOR_VERSION, [Minor version])
AC_DEFINE_UNQUOTED(MOP_MICRO_VERSION, $MOP_MICRO_VERSION, [Micro version])
AC_DEFINE_UNQUOTED(MOP_INTERFACE_AGE, $MOP_INTERFACE_AGE, [Interface age])
AC_DEFINE_UNQUOTED(MOP_BINARY_AGE, $MOP_BINARY_AGE, [Binary age])

dnl Initialize libtool
AM_PROG_LIBTOOL

dnl Initialize maintainer mode
AM_MAINTAINER_MODE

AC_PROG_CC

if test x$GCC = xyes ; then
  CFLAGS="$CFLAGS -Wall -Werror-implicit-function-declaration -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations"
fi

AC_PROG_AWK
AC_ISC_POSIX
AC_HEADER_STDC

AC_CHECK_PROG(GSL, gsl-config, yes, no)
if test x$GSL = xyes ; then
   gsl_cflags=`gsl-config --cflags`
   gsl_libs=`gsl-config --libs`
else
  echo \#################################################################
  echo \#
  echo \# GNU Scientific Library not available or not properly installed.
  echo \# You can get the most recent version from 
  echo \#      http://www.gnu.org/software/gsl/
  echo \#
  echo \#################################################################
  exit
fi

# taken from Stephane Popinet's GTS configure.in

AM_PATH_GLIB_2_0(2.4.0, has_glib2=yes, has_glib2=no, gthread gmodule)
if test x$has_glib2 = xyes; then
  glib_cflags=`$PKG_CONFIG glib-2.0 --cflags`
  glib_thread_cflags=`$PKG_CONFIG glib-2.0 --cflags gthread-2.0`
  glib_module_cflags=`$PKG_CONFIG glib-2.0 --cflags gmodule-2.0`
  glib_libs=`$PKG_CONFIG glib-2.0 --libs`
  glib_thread_libs=`$PKG_CONFIG glib-2.0 --libs gthread-2.0`
  glib_module_libs=`$PKG_CONFIG glib-2.0 --libs gmodule-2.0`
  GLIB_LIBS="$glib_libs"
  GLIB_DEPLIBS="$glib_libs"
else
  if test x$with_glib = xyes ; then
    AC_MSG_ERROR([
*** Directory must be specified for --with-glib])
  fi

  if test x$with_glib = x ; then 
    # Look for separately installed glib

    GLIB_REQUIRED_VERSION=1.2.8
    AM_PATH_GLIB($GLIB_REQUIRED_VERSION,,
      AC_MSG_ERROR([
*** GLIB $GLIB_REQUIRED_VERSION or better is required. The latest version of GLIB
*** is always available from ftp://ftp.gtk.org/.]),
      gmodule gthread)

    glib_cflags=`$GLIB_CONFIG --cflags`
    glib_thread_cflags=`$GLIB_CONFIG --cflags gthread`
    glib_module_cflags=`$GLIB_CONFIG --cflags gmodule`
    glib_libs=`$GLIB_CONFIG --libs`
    glib_thread_libs=`$GLIB_CONFIG --libs gthread`
    glib_module_libs=`$GLIB_CONFIG --libs gmodule`
    GLIB_LIBS="$glib_libs"
    GLIB_DEPLIBS="$glib_libs"
  else
    # Use uninstalled glib (assume they got the version right)

    GLIB_CONFIG=$with_glib/glib-config
    if test -x $GLIB_CONFIG ; then 
      :
    else
      AC_MSG_ERROR([GLIB directory ($with_glib) not present or not configured])
    fi

    # For use in mop-config
    glib_cflags=`$GLIB_CONFIG --cflags`
    glib_thread_cflags=`$GLIB_CONFIG --cflags gthread`
    glib_module_cflags=`$GLIB_CONFIG --cflags gmodule`
    glib_libs=`$GLIB_CONFIG --libs`
    glib_thread_libs=`$GLIB_CONFIG --libs gthread`
    glib_module_libs=`$GLIB_CONFIG --libs gmodule`

    glib_release=`$GLIB_CONFIG --version | sed 's%\\.[[0-9]]*$%%'`

    # canonicalize relative paths
    case $with_glib in 
      /*)
        glib_dir=$with_glib
        ;;
      *)
        glib_dir="\$(top_builddir)/$with_glib"
        ;;
    esac

    GLIB_CFLAGS="-I$glib_dir -I$glib_dir/gmodule"
    GLIB_LIBS="$glib_dir/libglib.la $glib_dir/gmodule.la"
    GLIB_DEPLIBS=

    AC_SUBST(GLIB_CFLAGS)
    AC_SUBST(GLIB_LIBS)
  fi
fi

AC_SUBST(glib_cflags)
AC_SUBST(glib_libs)
AC_SUBST(gsl_cflags)
AC_SUBST(gsl_libs)

CFLAGS="$CFLAGS $glib_cflags $gsl_cflags"
LDFLAGS="$LDFLAGS $glib_libs $gsl_libs"

AC_SUBST(CFLAGS)
AC_SUBST(CPPFLAGS)
AC_SUBST(LDFLAGS)

##AC_CHECK_LIB(m, sqrt)
AC_CONFIG_MACRO_DIR([m4])

LIBS="$glib_libs $gsl_libs $LIBS"

# Specify a configuration file
AC_CONFIG_HEADERS(config.h)

AC_CONFIG_FILES([
Makefile
tests/Makefile
doc/Makefile
doc/html/Makefile
src/Makefile
src/mop-config])

AC_OUTPUT
