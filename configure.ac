dnl tsshmm requires OpenMP and the general hidden Markov model (GHMM)
dnl library.
dnl
dnl GHMM does not yet distribute a pkg-config file, but the package includes
dnl it's own configuration tool ghmm-config:
dnl https://sourceforge.net/p/ghmm/code/HEAD/tree/trunk/ghmm/tools/ghmm-config.in
dnl This script will optionally use ghmm-config if present.  We cannot strictly
dnl rely on ghmm-config because a major operating system, Debian, neglects to
dnl include gmm-config tool with it's libghmm-dev package for some reason.

dnl Be careful not to introduce any trailing newlines in this tsshmm_VERSION m4
dnl macro.
m4_define(tsshmm_VERSION,
          m4_esyscmd([awk '/^Version/ {printf "%s", $NF; exit 0}' DESCRIPTION]))
AC_INIT(tsshmm, m4_defn([tsshmm_VERSION]))
AC_CONFIG_AUX_DIR([tools])
AC_CONFIG_MACRO_DIR([tools])
AM_INIT_AUTOMAKE
AM_PROG_AR
LT_INIT([disable-static])
# Use R-specific flags per
# https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Configure-and-cleanup
CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config --cppflags`
LDFLAGS=`"${R_HOME}/bin/R" CMD config LDFLAGS`
AC_LANG(C)
AC_OPENMP
AC_PROG_CC
AC_PROG_INSTALL

# Check for memory problems with an address sanitizer.
AC_ARG_ENABLE([asan],
              [AS_HELP_STRING([--enable-asan],
                              [When enabled, compile with -fsanitize=address.
                               Note that compilation and loading requires the LD_PRELOAD
                               environmental variable set, e.g.:
                               LD_PRELOAD=libasan.so.5
                               ...and for compilation and installation only,
                               disable ASAN with:
                               ASAN_OPTIONS=detect_leaks=0:use_odr_indicator=1
              ])] # AS_HELP_STRING
) # AC_ARG_ENABLE

AS_IF(test x"$enable_asan" = xyes,
      [AC_CHECK_LIB(asan, __sanitizer_report_error_summary,
                    CFLAGS="$CFLAGS -fsanitize=address"
                    LDFLAGS="$LDFLAGS -fsanitize=address"
      )]
)

# Read in any --with-ghmm user preference.
AC_ARG_WITH(ghmm,
            AS_HELP_STRING(--with-ghmm=GHMM_ROOT_OR_SEARCH_STRATEGY,
                           [path to GHMM installation, or a search strategy
                            of either "system", "bundled", or "fallback";
                            the last strategy first searches for GHMM in
                            system library paths -
                            including CPATH and LD_LIBRARY_PATH -
                            if that fails then fall back to using
                            the bundled GHMM library
                            @<:@default=fallback@:>@]
             )], # AS_HELP_STRING
             # action-if-given
             [],
             # action-if-not-given
             [: m4_divert_text([DEFAULTS], [with_ghmm=fallback])]
) # AC_ARG_WITH

# Overwrite the user preference --with-ghmm with the environmental variable
# GHMM_ROOT if set.
AC_MSG_CHECKING(for GHMM_ROOT environmental variable)
AS_IF(test "x$GHMM_ROOT" != x,
      dnl run-if-true
      with_ghmm=$GHMM_ROOT
      force_ghmm_root=yes
      AC_MSG_RESULT($GHMM_ROOT)
      AC_MSG_NOTICE(Using GHMM_ROOT ignores the --with-ghmm option),
      dnl run-if-false
      AC_MSG_RESULT(not set)
) dnl AS_IF

# Sanity check user input of --with-ghmm.
#
# Is --with-ghmm incorrectlly set to yes or no?
AS_IF(test x$force_ghmm_root != xyes,
      AC_MSG_CHECKING(whether --with-ghmm=GHMM_ROOT is incorrectly set to yes or no)
      [AS_CASE(x$with_ghmm,
      # Fail for these options
      xyes, [AC_MSG_FAILURE(
[Setting --with-ghmm or --with-ghmm=yes is insufficient.
Please specify "system", "bundled", "fallback" or the path to your GHMM
installation GHMM_ROOT either as a configure option:
    ./configure --with-ghmm=/path/to/your/ghmm/installation
or as an environmental variable:
    GHMM_ROOT=/path/to/your/ghmm/installation])],
      xno, [AC_MSG_FAILURE(
[Setting --without-ghmm or --with-ghmm=no is not allowed,
because GHMM is a required library.  Please specify "system", "bundled",
"fallback" or the path to your GHMM installation GHMM_ROOT either as a
configure option:
    ./configure --with-ghmm=/path/to/your/ghmm/installation
or as an environmental variable:
    GHMM_ROOT=/path/to/your/ghmm/installation])],
      # default
      AC_MSG_RESULT(neither yes nor no)
      )] # AS_CASE
) # AS_IF
#
# Is --with-ghmm set to a search strategy?
AS_IF(test x$force_ghmm_root != xyes,
      # run-if-true
      AC_MSG_CHECKING([for valid search strategy in --with-ghmm=STRATEGY])
      [AS_CASE(x$with_ghmm,
               # Allow these options
               xsystem, with_ghmm_strategy=system,
               xbundled, with_ghmm_strategy=bundled,
               xfallback, with_ghmm_strategy=fallback,
               # default
               with_ghmm_strategy=no
      )] # AS_CASE
      [AS_IF(test x$with_ghmm_strategy != xno,
             AC_MSG_RESULT($with_ghmm_strategy),
             [AC_MSG_RESULT([invalid, but might be --with-ghmm=GHMM_ROOT])]
      )] # AS_IF
) # AS_IF
#
# Is GHMM_ROOT from the environmental variable or from --with-ghmm a valid
# path?  Either is a valid path invalid path or invalid input.
AS_IF([test x$with_ghmm_strategy = x],
      AC_MSG_NOTICE(checking whether GHMM_ROOT=$with_ghmm is a valid path)
      AS_IF([test -d $with_ghmm],
            [with_ghmm_strategy=no],
            AC_MSG_FAILURE(The directory GHMM_ROOT=$with_ghmm was not found.)
      ) # AS_IF
) # AS_IF

AC_SUBST(GHMM_LIBS, -lghmm)

# If GHMM_ROOT was provided, set the header and library paths.
#
# Check for the existance of include/ and lib/ sub-directories and if both are
# present add them to CFLAGS and LDFLAGS per
# https://stackoverflow.com/a/44534707
AS_IF(test x$with_ghmm_strategy = xno,
      AS_IF([test -f $with_ghmm/include/ghmm/ghmm.h],
            [AC_SUBST(GHMM_CPPFLAGS, ["-I$with_ghmm/include"])]
            AC_MSG_FAILURE(GHMM header not found.))
      AS_IF([test -f $with_ghmm/lib/libghmm.so],
            [AC_SUBST(GHMM_LDFLAGS,
                     ["-Wl,-rpath=$with_ghmm/lib -L$with_ghmm/lib"])],
            AC_MSG_FAILURE(GHMM library not found.))
) # AS_IF

AS_IF(test x$with_ghmm_strategy = xsystem -o \
           x$with_ghmm_strategy = xfallback,
      AS_IF(test x$with_ghmm_strategy = xsystem,
            AC_MSG_NOTICE(Checking GHMM only on system),
            AC_MSG_NOTICE(Checking GHMM on system))
      # Check for system or LD_LIBRARY_PATH installed GHMM.
      AC_CHECK_LIB(ghmm, ghmm_dmodel_check, found_ghmm_system_lib=yes)
      # Check for system or CPATH / C_INCLUDE_PATH GHMM.
      [AC_CHECK_HEADERS(ghmm/ghmm.h, found_ghmm_system_headers=yes)]
      AS_IF(test x$found_ghmm_system_lib = xyes -a \
                 x$found_ghmm_system_headers = xyes,
            found_ghmm_system=yes
            AC_MSG_NOTICE(GHMM found on system.),
            AS_IF(test x$with_ghmm_strategy = xsystem,
                  AC_MSG_FAILURE(GHMM library not found.)))
) # AS_IF

# Use bundled dependency if necessary.
AM_CONDITIONAL(BUNDLED_GHMM, false)
AS_IF(test x$found_ghmm_system != xyes &&
      test x$with_ghmm_strategy = xbundled -o \
           x$with_ghmm_strategy = xfallback,
      AS_IF(test x$with_ghmm_strategy = xbundled,
            AC_MSG_NOTICE(Using bundled GHMM),
            AC_MSG_NOTICE(Falling back to using bundled GHMM))
      AM_CONDITIONAL(BUNDLED_GHMM, true)
      [AX_SUBDIRS_CONFIGURE([src/ghmm-0.9-rc3],
                            [[CFLAGS=$CFLAGS],
			     [--enable-gsl=no],
                             [--disable-gsltest],
                             [--with-rng=mt],
                             [--with-python=no],
                             [--libdir=$PWD/libs],
                             [--prefix=$PWD/src/libunix]])]
      # Our recursive automake rule needs the absolute path of GHMM_ROOT.
      [AC_SUBST(GHMM_ROOT, ["$PWD/src/ghmm-0.9-rc3/ghmm"])]
      [AS_IF([test -d $GHMM_ROOT], [],
             AC_MSG_FAILURE(Directory of bundled GHMM not found.))]
      [AC_SUBST(GHMM_CPPFLAGS, ["-I$GHMM_ROOT/.."])]
      # Using -rpath=. prefers the bundled over any system installation.
      [AC_SUBST(GHMM_LDFLAGS, ["-Wl,-rpath=. -L$GHMM_ROOT/.libs"])]
      AC_MSG_NOTICE(Applying patches to GHMM to fix errors and warnings from "R CMD check")
      # Patch bug in upstream's configure bug:
      #
      #     config.status: executing libtool commands
      #     /usr/bin/rm: cannot remove 'libtoolT': No such file or directory
      #
      # ... it's not yet clear how to fix upstream's configure.ac.  The command
      # below is idempotent due to the regex address match prefix, so rerunning
      # the patch when rerunning ./configure has no adverse effect.
      [cp -au src/ghmm-0.9-rc3/configure{,.bak}]
      sed -i "/RM='@<:@^f@:>@*'/s/'$/ -f'/" src/ghmm-0.9-rc3/configure
      [touch -r src/ghmm-0.9-rc3/configure{.bak,}]
      [diff -u src/ghmm-0.9-rc3/configure{.bak,}]
      # Patch bug in upstream's xmlwriter.c when compiling with
      # -Werror=format-security:
      #
      #     xmlwriter.c: In function 'writeDiscreteStateContents':
      #     xmlwriter.c:356:3: error: format not a string literal and no format arguments [-Werror=format-security]
      #        rc = xmlTextWriterWriteFormatElement(writer, "backgroundKey", f->model.d[moNo]->bp->name[bgId]);
      #        ^~The command
      #
      # This command is also idempotent, therefore rerunning ./configure has no
      # adverse effect.
      [cp -au src/ghmm-0.9-rc3/ghmm/xmlwriter.c{,.bak}]
      [sed -i '356s/y", f/y", "%s", f/' src/ghmm-0.9-rc3/ghmm/xmlwriter.c]
      [touch -r src/ghmm-0.9-rc3/ghmm/xmlwriter.c{.bak,}]
      [diff -u src/ghmm-0.9-rc3/ghmm/xmlwriter.c{.bak,}]
      # Patch bug in upstreams mcmc.c test when compiling with -Wint-conversion:
      #
      #     mcmc.c:146:15: warning: implicit declaration of function ‘ghmm_dmodel_cfbgibbs’; did you mean ‘ghmm_dmodel_logp’? [-Wimplicit-function-declaration]
      #        int * Q =   ghmm_dmodel_cfbgibbs(mo, my_output,
      #                    ^~~~~~~~~~~~~~~~~~~~
      #                    ghmm_dmodel_logp
      #     mcmc.c:146:15: warning: initialization of ‘int *’ from ‘int’ makes pointer from integer without a cast [-Wint-conversion]
      #     mcmc.c:148:108: warning: passing argument 4 of ‘ghmm_dmodel_viterbi_logp’ makes pointer from integer without a cast [-Wint-conversion]
      #        printf("viterbi prob mcmc%f \n", ghmm_dmodel_viterbi_logp(mo, my_output->seq[0], my_output->seq_len[0], Q[0]));
      #
      # This command is also idempotent, therefore rerunning ./configure has no
      # adverse effect.
      [cp -au src/ghmm-0.9-rc3/tests/mcmc.c{,.bak}]
      sed -i 's/int @<:@*@:>@ Q/int **Q/' src/ghmm-0.9-rc3/tests/mcmc.c
      grep -q 'ghmm@<:@/@:>@cfbgibbs@<:@.@:>@h' src/ghmm-0.9-rc3/tests/mcmc.c || sed -i '9i\
#include <ghmm/cfbgibbs.h>\
#include <ghmm/fbgibbs.h>' src/ghmm-0.9-rc3/tests/mcmc.c
      [touch -r src/ghmm-0.9-rc3/tests/mcmc.c{.bak,}]
      [diff -u src/ghmm-0.9-rc3/tests/mcmc.c{.bak,}]
      AC_MSG_NOTICE(Finished patching GHMM)
) # AS_IF

# Variables for Doxygen.
AC_SUBST(R_INCLUDE, ["${CPPFLAGS#-I}"])
AC_SUBST(GHMM_INCLUDE, ["${GHMM_CPPFLAGS#-I}"])

AS_ECHO
AS_BOX(Configuration of $PACKAGE_STRING)
AS_ECHO
AS_ECHO("             CC = $CC")
AS_ECHO
AS_ECHO("         CFLAGS = $CFLAGS")
AS_ECHO("       CPPFLAGS = $CPPFLAGS")
AS_ECHO("        LDFLAGS = $LDFLAGS")
AS_ECHO("           LIBS = $LIBS")
AS_ECHO
AS_ECHO("  OPENMP_CFLAGS = $OPENMP_CFLAGS")
AS_ECHO
AS_ECHO("      GHMM_ROOT = $GHMM_ROOT")
AS_ECHO("  GHMM_CPPFLAGS = $GHMM_CPPFLAGS")
AS_ECHO("   GHMM_LDFLAGS = $GHMM_LDFLAGS")
AS_ECHO("      GHMM_LIBS = $GHMM_LIBS")
AS_ECHO
AS_ECHO("          CPATH = $CPATH")
AS_ECHO(" C_INCLUDE_PATH = $C_INCLUDE_PATH")
AS_ECHO("LD_LIBRARY_PATH = $LD_LIBRARY_PATH")
AS_ECHO
AS_ECHO("      R_INCLUDE = $R_INCLUDE")
AS_ECHO("   GHMM_INCLUDE = $GHMM_INCLUDE")
AS_ECHO

AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([src/Doxyfile])
AC_OUTPUT
