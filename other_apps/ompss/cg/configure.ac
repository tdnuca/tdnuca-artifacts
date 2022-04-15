#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([ompss_cg], [1.1], [luc.jaulmes@bsc.es])
AC_OUTPUT(Makefile src/Makefile)
AC_CONFIG_SRCDIR([src])
AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([no-dependencies])
AC_LANG([C])


# now more custom stuff:
# get _HOME from environment for all (potential) prerequesites)
AC_ARG_VAR([EXTRAE_HOME], [Extrae installation path])
AC_ARG_VAR([NANOS_HOME], [Nanos installation path])
AC_ARG_VAR([MCC_HOME], [Mercurium installation path])
AC_ARG_VAR([MCC], [Mercurium C compiler command (default=mcc)])
AC_ARG_VAR([MPI_HOME], [MPI installation path])


# allow overrides with --package{,-lib,-include}=DIR
# and define PACKAGE_LIBDIR=/path/to/libs PACKAGE_INCLUDEDIR=/path/to/headers

# extrae package, has a default of no
AC_ARG_WITH([extrae],
	[AS_HELP_STRING([--with-extrae[[=DIR]]], [use extrae to emit events into traces rather than printing to stdout, optionally specify path])],
	[AS_IF([test "x$withval" != "xyes"],
		[EXTRAE_HOME=$withval]
	)],
	[with_extrae="no"] # if flag not given, disable extrae
)

AC_ARG_WITH([extrae-lib],
	[AS_HELP_STRING([--with-extrae-lib=DIR], [override path to extrae libraries (only with --with-extrae)])],
	[EXTRAE_LIBDIR=$withval],
	[AS_IF([test -d "$EXTRAE_HOME"],
		[AS_IF([test -d "$EXTRAE_HOME/lib64"],
			[EXTRAE_LIBDIR=$EXTRAE_HOME/lib64],
			[EXTRAE_LIBDIR=$EXTRAE_HOME/lib]
		)],
		[AC_MSG_NOTICE([EXTRAE_HOME=$EXTRAE_HOME but -n false])]
	)]
)

AC_ARG_WITH([extrae-include],
	[AS_HELP_STRING([--with-extrae-include=DIR], [override path to extrae headers (only with --with-extrae)])],
	[EXTRAE_INCLUDEDIR=$withval],
	[AS_IF([test -d "$EXTRAE_HOME"],
		[EXTRAE_INCLUDEDIR=$EXTRAE_HOME/include]
	)]
)


# nanos
AC_ARG_WITH([nanos],
	[AS_HELP_STRING([--with-nanos=DIR], [override path to nanos installation])],
	[AC_SUBST([NANOS_HOME], [$withval])]
)

AC_ARG_WITH([nanos-lib],
	[AS_HELP_STRING([--with-nanos-lib=DIR], [override path to nanos libraries])],
	[LDFLAGS+=" -L$withval"],
	[AS_IF([test -d "$NANOS_HOME/lib64"],
		[LDFLAGS+=" -L$NANOS_HOME/lib64"],
		[AS_IF([test -d "$NANOS_HOME/lib"],
			[LDFLAGS+=" -L$NANOS_HOME/lib"]
		)]
	)]
)

AC_ARG_WITH([nanos-include],
	[AS_HELP_STRING([--with-nanos-include=DIR], [override path to nanos headers])],
	[CPPFLAGS+=" -I$withval"],
	[AS_IF([test -d "$NANOS_HOME"],
		[CPPFLAGS+=" -I$NANOS_HOME/include"]
	)]
)

# MPI
AC_ARG_WITH([mpi],
	[AS_HELP_STRING([--with-mpi=DIR], [override path to mpi installation])],
	[AS_IF([test -d "$withval"],
		[AC_SUBST([MPI_HOME], [$withval])]
	)]
)

AC_ARG_WITH([mpi-lib],
	[AS_HELP_STRING([--with-mpi-lib=DIR], [override path to mpi libraries])],
	[LDFLAGS+=" -L$withval"],
	[AS_IF([test -d "$MPI_HOME/lib64"],
		[LDFLAGS+=" -L$MPI_HOME/lib64"],
		[AS_IF([test -d "$MPI_HOME/lib"],
			[LDFLAGS+=" -L$MPI_HOME/lib"]
		)]
	)]
)

AC_ARG_WITH([mpi-include],
	[AS_HELP_STRING([--with-mpi-include=DIR], [override path to mpi headers])],
	[CPPFLAGS+=" -I$withval"],
	[AS_IF([test -d "$MPI_HOME"],
		[CPPFLAGS+=" -I$MPI_HOME/include"]
	)]
)

# mercurium
AC_ARG_WITH([mercurium],
	[AS_HELP_STRING([--with-mercurium=DIR], [override path to the mercurium compiler])],
	[AC_SUBST([MCC_HOME], [$withval])]
)

AS_IF([test "x$with_mpi" == xno], [
   AC_MSG_NOTICE([Building OmpSs-only CG])
], [
   AC_MSG_NOTICE([Building MPI+OmpSs CG])
   AC_DEFINE([USE_MPI], [1], [Whether to build a hybrid MPI+OmpSs or single-node OmpSs CG])
])

# Don't let user say --without-{nanos,mercurium}
AS_IF([test "x$with_nanos" == xno], [AC_MSG_ERROR([Nanos is required to build CG])])
AS_IF([test "x$with_mercurium" == xno], [AC_MSG_ERROR([Mercurium is required to build CG])])

# if --with-extrae was given (not --without-extrae nor --with-extrae=no) then define USE_EXTRAE
AM_CONDITIONAL([USE_EXTRAE], [test "x$with_extrae" != xno])


# Find a backend compiler.
AC_PROG_CC([gcc icc cc])

# Checks for libraries.
AC_CHECK_LIB(m, sqrt)

# Checks for header files.
AC_CHECK_HEADERS([float.h stdlib.h string.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE

# Checks for library functions.
AC_HEADER_MAJOR
AC_FUNC_MALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([gettimeofday strtol])
AC_PROG_RANLIB

AC_CHECK_HEADER([nanox/nanos_omp.h], [],
	[AC_MSG_ERROR([Cannot find nanos headers, please specify --with-nanos or --with-nanos-includes])]
)

AS_IF([test "x$with_mpi" != xno], [
	AC_CHECK_HEADER([mpi.h], [],
		[AC_MSG_ERROR([Cannot find MPI headers, please specify --with-mpi or --with-mpi-includes])]
	)
	AC_SEARCH_LIBS([MPI_Init], [mpi mpich], [], [AC_MSG_FAILURE([Cannot find MPI_Init in any library])])
])

# Check that mcc exists, with command given through MCC. If found, set in CC_FOR_OMPSS
# and check its functionality by setting *_FOR_OMPSS and running AC_PROG_CC_FOR_OMPSS
AS_IF([test -z "$MCC"],
	[MCC=mcc]
)
AC_MSG_NOTICE([Using MCC=$MCC for the mercurium compiler])
AS_IF([test -z "$MCC_HOME"],
	[AC_CHECK_PROG([CC_FOR_OMPSS], [$MCC], [$MCC], [not found])],
	[AC_PATH_PROG([CC_FOR_OMPSS], [$MCC], [not found], [$MCC_HOME/bin$PATH_SEPARATOR$PATH])]
)
AS_IF([test "x$CC_FOR_OMPSS" == "xnot found"],
	[AC_MSG_ERROR([mercurium compiler $MCC cannot be found, please install it or specify its name via MCC, and/or its path via MCC_HOME or --with-mercurium=DIR])],
)

CPP_FOR_OMPSS="$CC_FOR_OMPSS -E"
CFLAGS_FOR_OMPSS="--cc=$CC $CFLAGS"
CPPFLAGS_FOR_OMPSS="--cpp=$CC $CPPFLAGS"
LDFLAGS_FOR_OMPSS="--ld=$CC $LDFLAGS"

AC_PROG_CC_FOR_OMPSS

AS_IF([test "x$with_mpi" == xno],
	[seqtrace="seqtrace"],
	[seqtrace="mpitrace"]
)
AS_IF([test "x$with_mpi" == xno],
	[nanostrace="nanostrace"],
	[nanostrace="nanosmpitrace"]
)

# If extrae is enabled, we need to check linking libnanostrace and libseqtrace
AS_IF([test "x$with_extrae" != xno], [
	AC_DEFINE([EXTRAE_EVENTS], [1], [Define to use custom extrae events])
	AS_IF([test -d "$EXTRAE_INCLUDEDIR"], [CPPFLAGS+=" -I$EXTRAE_INCLUDEDIR"] [CPPFLAGS_FOR_OMPSS+=" -I$EXTRAE_INCLUDEDIR"])
	AS_IF([test -d "$EXTRAE_LIBDIR"],     [LDFLAGS+=" -L$EXTRAE_LIBDIR"] [LDFLAGS_FOR_OMPSS+=" -L$EXTRAE_LIBDIR"])
	AC_CHECK_LIB(["$seqtrace"], [main],
		[
			AC_SUBST([LIBEXTRAE_SEQ], ["-l$seqtrace"])
		],
		[AS_IF([test "x$with_extrae" != xcheck],
			[AC_MSG_FAILURE([--with-extrae was given but linking with lib$seqtrace failed])]
		)]
	)
	AC_CHECK_LIB(["$nanostrace"], [main],
		[
			AC_SUBST([LIBEXTRAE_NANOS], ["-l$nanostrace"])
		],
		[AS_IF([test "x$with_extrae" != xcheck],
			[AC_MSG_FAILURE([--with-extrae was given but linking with lib$nanostrace failed])]
		)]
	)
])

# Set verbosity level
AC_ARG_ENABLE([verbosity], [AS_HELP_STRING([--enable-verbosity[[=n]]], [enable verbosity of CG and set level (default is disabled)])])

AS_CASE(["x$enable_verbosity"],
	[x[[1-4]]], [AC_DEFINE([VERBOSE], $enable_verbosity, [Set level of verbosity])],
	[xyes], [AC_DEFINE([VERBOSE], [1], [Set level of verbosity])],
	[AC_DEFINE([PERFORMANCE], [1], [Set non-verbose])]
)

AC_ARG_ENABLE([reductions], [AS_HELP_STRING([--enable-reductions[[=n]]], [enable task-reductions in CG, instead of using concurrent + atomic (see README)])])
AS_CASE(["x$enable_reductions"],
	[xyes], [AC_DEFINE([USE_REDUCTIONS], [1], [Use task reductions])],
	[xno], [],
	[AC_MSG_ERROR([Must select to enable or disable reductions.])]
)

# Output compiler as backend, and mcc as compiler
AC_SUBST([CC_BACKEND], [$CC])
AC_SUBST([CFLAGS_BACKEND], [$CFLAGS])
AC_SUBST([CPPFLAGS_BACKEND], [$CPPFLAGS])
AC_SUBST([LDFLAGS_BACKEND], [$LDFLAGS])

AC_SUBST([CC], [$CC_FOR_OMPSS])
AC_SUBST([CFLAGS], [$CFLAGS_FOR_OMPSS])
AC_SUBST([CPPFLAGS], [$CPPFLAGS_FOR_OMPSS])
AC_SUBST([LDFLAGS], [$LDFLAGS_FOR_OMPSS])

AC_OUTPUT

