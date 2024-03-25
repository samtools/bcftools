# Configure script for bcftools, utilities for VCF/BCF files.
#
#    Copyright (C) 2015,2017-2020 Genome Research Ltd.
#
#    Author: John Marshall <jm18@sanger.ac.uk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

dnl Process this file with autoconf to produce a configure script
AC_INIT([BCFtools], m4_esyscmd_s([./version.sh 2>/dev/null]),
        [samtools-help@lists.sourceforge.net], [], [http://www.htslib.org/])
AC_PREREQ([2.63])  dnl This version introduced 4-argument AC_CHECK_HEADER
AC_CONFIG_SRCDIR([main.c])
AC_CONFIG_HEADERS([config.h])

m4_include([m4/ax_with_htslib.m4])
m4_include([m4/hts_prog_cc_warnings.m4])

dnl Copyright notice to be copied into the generated configure script
AC_COPYRIGHT([Portions copyright (C) 2015,2017-2020 Genome Research Ltd.

This configure script is free software: you are free to change and
redistribute it.  There is NO WARRANTY, to the extent permitted by law.])

AC_PROG_CC

dnl Turn on compiler warnings, if possible
HTS_PROG_CC_WARNINGS
dnl Flags to treat warnings as errors.  These need to be applied to CFLAGS
dnl later as they can interfere with some of the tests (notably AC_SEARCH_LIBS)
HTS_PROG_CC_WERROR(hts_late_cflags)

AC_SYS_LARGEFILE

AC_ARG_ENABLE([bcftools-plugins],
  [AS_HELP_STRING([--disable-bcftools-plugins],
                  [do not build optional plugins])],
  [], [enable_bcftools_plugins=yes])
AC_SUBST(enable_bcftools_plugins)

AC_ARG_WITH([bcf-plugin-dir],
  [AS_HELP_STRING([--with-bcf-plugin-dir=DIR],
                  [plugin installation location [LIBEXECDIR/bcftools]])],
  [case $withval in
     yes|no) AC_MSG_ERROR([no directory specified for --with-plugin-dir]) ;;
   esac],
   [with_bcf_plugin_dir='$(libexecdir)/bcftools'])
AC_SUBST([bcf_plugindir], $with_bcf_plugin_dir)

AC_ARG_WITH([bcf-plugin-path],
  [AS_HELP_STRING([--with-bcf-plugin-path=PATH],
                  [default BCFTOOLS_PLUGINS plugin search path [PLUGINDIR]])],
  [case $withval in
     yes) AC_MSG_ERROR([no path specified for --with-bcf-plugin-path]) ;;
     no)  with_bcf_plugin_path= ;;
   esac],
  [with_bcf_plugin_path=$with_bcf_plugin_dir])
AC_SUBST([bcf_pluginpath], $with_bcf_plugin_path)

AC_ARG_ENABLE([perl-filters],
  [AS_HELP_STRING([--enable-perl-filters],
                  [build in support for PERL scripts in -i/-e filtering expressions])],
  [], [enable_perl_filters=no])
AC_SUBST(enable_perl_filters)

AC_ARG_ENABLE([libgsl],
  [AS_HELP_STRING([--enable-libgsl],
                  [build options that require the GNU scientific library (changes bcftools license to GPL, see INSTALL for details)])],
  [], [enable_libgsl=no])

AC_ARG_WITH([cblas],
  [AS_HELP_STRING([--with-cblas],
                  [which cblas library to use (only needed with libgsl; default: try cblas then gslcblas)])],
  [AS_CASE([$withval],
	   [no],[AC_MSG_ERROR([no library name specified for --with-cblas])])],
  [with_cblas=yes])

AX_WITH_HTSLIB
if test "$ax_cv_htslib" != yes; then
  AC_MSG_ERROR([HTSlib development files not found

BCFtools uses HTSlib to parse bioinformatics file formats etc.  Building it
requires an unpacked HTSlib source tree (which will be built in conjunction
with bcftools) or a previously-installed HTSlib.  In either case you may
need to configure --with-htslib=DIR to locate the appropriate HTSlib.

FAILED.  You must supply an HTSlib in order to build bcftools successfully.])
fi

if test "$ax_cv_htslib_which" = source; then
  Hsource=
  Hinstall='#'
else
  Hsource='#'
  Hinstall=
fi
AC_SUBST([Hsource])
AC_SUBST([Hinstall])

AS_IF([test "$enable_bcftools_plugins" != "no"], [dnl
  AC_MSG_CHECKING([shared library type])
  test -n "$host_alias" || host_alias=unknown-`uname -s`
  AS_CASE([$host_alias],
    [*-cygwin* | *-CYGWIN*],[dnl
      host_result="Cygwin DLL"
      PLATFORM=CYGWIN
      PLUGIN_EXT=.cygdll
      HTSLIB_DLL=libhts.dll.a
	  ],

    [*-msys* | *-MSYS* | *-mingw* | *-MINGW*],[dnl
      host_result="MSYS dll"
      PLATFORM=MSYS
      PLUGIN_EXT=.dll
      HTSLIB_DLL=hts.dll.a
      # This also sets __USE_MINGW_ANSI_STDIO which in turn makes PRId64,
      # %lld and %z printf formats work.  It also enforces the snprintf to
      # be C99 compliant so it returns the correct values (in kstring.c).
      CPPFLAGS="$CPPFLAGS -D_XOPEN_SOURCE=600"],

    [*-darwin* | *-Darwin*],[dnl
      host_result="Darwin dylib"
      PLATFORM=Darwin
      PLUGIN_EXT=.so],

    [*],[dnl
      host_result="plain .so"
      PLATFORM=default
      PLUGIN_EXT=.so])
  AC_MSG_RESULT([$host_result])
  AC_SUBST([PLATFORM])
  AC_SUBST([PLUGIN_EXT])
  AC_DEFINE_UNQUOTED([PLUGIN_EXT], ["$PLUGIN_EXT"],
                     [Platform-dependent plugin filename extension.])
  AC_SUBST([HTSLIB_DLL])
  AC_DEFINE([ENABLE_BCF_PLUGINS], 1,
            [Define if BCFtools should enable plugins.])

  AS_CASE([$PLUGIN_EXT],
    [*cygdll | *so],
	  [AC_SEARCH_LIBS([dlopen], [dl], [],
        [AC_MSG_ERROR([dlopen() not found
Plugin support requires dynamic linking facilities from the operating system.
Either configure with --disable-bcftools-plugins or resolve this error to
build bcftools.])])],
    [*dll],
	  [AC_SEARCH_LIBS([LoadLibraryA], [], [],
        [AC_MSG_ERROR([LoadLibraryA() not found
Plugin support requires dynamic linking facilities from the operating system.
Either configure with --disable-bcftools-plugins or resolve this error to
build bcftools.])])],
	)

  AC_MSG_CHECKING([if the compiler accepts -rdynamic])
  save_CFLAGS=$CFLAGS
  CFLAGS="$CFLAGS -rdynamic"
  AC_LINK_IFELSE([AC_LANG_PROGRAM()],[cc_rdynamic=yes],[cc_rdynamic=no])
  CFLAGS=$save_CFLAGS
  AC_MSG_RESULT([$cc_rdynamic])
  AS_IF([test "$cc_rdynamic" = "yes"],dnl
    [CC_RDYNAMIC_OPT=-rdynamic],[CC_RDYNAMIC_OPT=])
  AC_SUBST([CC_RDYNAMIC_OPT])
])
save_LIBS=$LIBS
zlib_devel=ok
dnl Set a trivial non-empty INCLUDES to avoid excess default includes tests
AC_CHECK_HEADER([zlib.h], [], [zlib_devel=missing], [;])
AC_CHECK_LIB(z, inflate,  [], [zlib_devel=missing])
LIBS=$save_LIBS

if test $zlib_devel != ok; then
  AC_MSG_ERROR([zlib development files not found

BCFtools uses compression routines from the zlib library <http://zlib.net>.
Building bcftools requires zlib development files to be installed on the build
machine; you may need to ensure a package such as zlib1g-dev (on Debian or
Ubuntu Linux) or zlib-devel (on RPM-based Linux distributions) is installed.

FAILED.  This error must be resolved in order to build bcftools successfully.])
fi

AS_IF([test "$enable_libgsl" != "no"], [dnl
  gsl_devel=ok
  save_LIBS="$LIBS"dnl Libraries for GSL are listed in $(GSL_LIBS)
  GSL_LIBS=
  AC_CHECK_HEADER([gsl/gsl_version.h], [], [gsl_devel=missing], [;])
  AS_IF([test "$gsl_devel" != "missing"], [dnl
    AS_CASE([$with_cblas],
      [yes],[AC_SEARCH_LIBS([cblas_dgemm],[cblas gslcblas],
              [GSL_LIBS="$ac_cv_search_cblas_dgemm"],
              [gsl_devel=missing])],
      [AC_CHECK_LIB([$with_cblas],[cblas_dgemm],
         [GSL_LIBS="-l$with_cblas"],[gsl_devel=missing])])
  ])
  AS_IF([test "$gsl_devel" != "missing"], [dnl
    AC_CHECK_LIB(gsl, gsl_multifit_gradient,
      [GSL_LIBS="-lgsl $GSL_LIBS"], [gsl_devel=missing], [$GSL_LIBS])
  ])

  AS_IF([test "$gsl_devel" != "ok"], [dnl
    AC_MSG_ERROR([libgsl development files not found

The bcftools polysomy command uses routines from the GNU Scientific library
<https://www.gnu.org/software/gsl/>.  Building BCFtools with libgsl enabled
requires libgsl development files to be installed on the build machine; you
may need to ensure a package such as libgsl0-dev (on Debian or Ubuntu linux)
or gsl-devel (on RPM-based Linux distributions) is installed.

Either configure with --disable-libgsl or resolve this error to build bcftools.
    ])
  ])

  LIBS="$save_LIBS"
  USE_GPL=1
  AC_SUBST([USE_GPL])
  AC_SUBST([GSL_LIBS])
])

AS_IF([test "$enable_perl_filters" != "no" ], [dnl
    save_LIBS="$LIBS"
    if ! perl -MExtUtils::Embed -e ccopts -e ldopts >/dev/null 2>&1; then
       AC_MSG_ERROR([Cannot `perl -MExtUtils::Embed -e ccopts -e ldopts`])
    fi
    AC_DEFINE([ENABLE_PERL_FILTERS], 1, [Define if BCFtools should enable for support PERL scripts in -i/-e filtering expressions.])
    PERL_CCOPTS="`perl -MExtUtils::Embed -e ccopts`"
    PERL_LIBS="`perl -MExtUtils::Embed -e ldopts`"
    dnl A hack to remove faulty "-arch i386" on mac as it breaks linking
    if test "x$PLATFORM" = "xDarwin"; then
        PERL_CCOPTS="`perl -MExtUtils::Embed -e ccopts | sed 's,-arch i386,,'`"
        PERL_LIBS="`perl -MExtUtils::Embed -e ldopts | sed 's,-arch i386,,'`"
    fi
    AS_CASE([$PERL_LIBS],
      [*redhat/redhat-hardened-ld*],[[
        PERL_CCOPTS="`echo $PERL_CCOPTS | sed 's,-specs=[a-z/]*/redhat/redhat-hardened-[a-z0-9]* *,,g'`"
        PERL_LIBS="`echo $PERL_LIBS | sed 's,-specs=[a-z/]*/redhat/redhat-hardened-[a-z0-9]* *,,g'`"]
        AC_MSG_NOTICE([removing -specs=...redhat-hardened... from Perl options])])
    AC_SUBST([PERL_CCOPTS])
    AC_SUBST([PERL_LIBS])

    dnl Look for libperl
    AC_SEARCH_LIBS(Perl_call_pv, perl, [], [
        AC_MSG_ERROR([Unable to find the Perl_call_pv() function for perl filters])
    ],[$PERL_CCOPTS $PERL_LIBS])
    LIBS="$save_LIBS"
])

dnl Apply value from HTS_PROG_CC_WERROR (if set)
AS_IF([test "x$hts_late_cflags" != x],[CFLAGS="$CFLAGS $hts_late_cflags"])

dnl Look for regcomp in various libraries (needed on windows/mingw).
AC_SEARCH_LIBS(regcomp, regex, [libregex=needed], [])

dnl Force POSIX mode on Windows/Mingw
test -n "$host_alias" || host_alias=unknown-`uname -s`
case $host_alias in
  *-msys* | *-MSYS* | *-mingw* | *-MINGW*)
    host_result="MSYS dll"
    PLATFORM=MSYS
    PLUGIN_EXT=.dll
    # This also sets __USE_MINGW_ANSI_STDIO which in turn makes PRId64,
    # %lld and %z printf formats work.  It also enforces the snprintf to
    # be C99 compliant so it returns the correct values (in kstring.c).
    CPPFLAGS="$CPPFLAGS -D_XOPEN_SOURCE=600"
    ;;
esac

AC_CONFIG_FILES([config.mk])
AC_OUTPUT

AS_IF([test "x$USE_GPL" != "x"], [dnl
  AC_MSG_NOTICE([License has been changed to GPL.

Linking against the GNU scientific library requires that the resulting
binary is licensed using the GNU General Public License.  In particular,
this means you will have to follow the terms of this license in order
to distribute binary copies.  For more information see the INSTALL file.

Text of the GNU General Public License can be found in the LICENSE file,
or at http://www.gnu.org/copyleft/gpl.html

Information on the GNU Scientific Library, including the rationale for
its choice of license can be found at https://www.gnu.org/software/gsl/
])])
