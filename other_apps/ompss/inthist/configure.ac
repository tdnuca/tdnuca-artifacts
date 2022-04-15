#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.64)
AC_INIT([inthist_crossweave],[1.1.1],[isaac.sanchez@bsc.es (orig: pieter.bellens@bsc.es)])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/inthist_main.c])
AM_INIT_AUTOMAKE
AC_CONFIG_HEADER([src/config.h])
AC_CONFIG_FILES([Makefile 
		src/Makefile])
AC_CONFIG_MACRO_DIR([m4])
AC_LANG([C])

: ${CFLAGS="-Wall"}
: ${CFLAGS_PERF="-O3"}
: ${CFLAGS_INSTR="-O3"}
: ${CFLAGS_DEBUG="-g"}
: ${LDFLAGS_PERF="-O3"}
: ${LDFLAGS_INSTR=""}
: ${LDFLAGS_DEBUG="-g"}
: ${LIBS="-lm"}

AC_ARG_ENABLE([ompss],
	AS_HELP_STRING([--disable-ompss], [disable OmpSs; enabled by default]))

case "x$enable_ompss" in
xno)
	AC_DEFINE([USE_OMPSS], 0, [Disable OmpSs])
	test_CC="gcc icc xlc cc"
;;
*)
	AC_DEFINE([USE_OMPSS], 1, [Enable OmpSs])
	enable_ompss="yes"
	test_CC="smpcc mcc imcc xlmcc"
esac

AC_PROG_CC([$test_CC])

# Checks for libraries.
AM_PROG_CC_C_O

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h malloc.h stdlib.h string.h sys/time.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T
AC_HEADER_TIME

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([floor gettimeofday memset strchr strstr])

case "x$enable_ompss" in
xyes)
	AX_CHECK_COMPILE_FLAG([--ompss],
		[AX_APPEND_FLAG([--ompss], [CFLAGS])],
		[AC_MSG_ERROR([--ompss not supported. Are you sure this is an OmpSs compiler?])])

	AX_CHECK_COMPILE_FLAG([--ompss],
		[AX_APPEND_FLAG([--ompss], [LDFLAGS])],
		[AC_MSG_ERROR([--ompss not supported. Are you sure this is an OmpSs compiler?])])

	AX_CHECK_COMPILE_FLAG([--instrument],
		[AX_APPEND_FLAG([--instrument], [CFLAGS_INSTR])],
		[AC_MSG_WARN([--instrument not supported. Are you sure this is an OmpSs compiler?])])

	AX_CHECK_COMPILE_FLAG([--instrument],
		[AX_APPEND_FLAG([--instrument], [LDFLAGS_INSTR])],
		[AC_MSG_WARN([--instrument not supported. Are you sure this is an OmpSs compiler?])])

	AX_CHECK_COMPILE_FLAG([--debug],
		[AX_APPEND_FLAG([--debug -k], [CFLAGS_DEBUG])],
		[AC_MSG_WARN([--debug not supported. Are you sure this is an OmpSs compiler?])])

	AX_CHECK_COMPILE_FLAG([--debug],
		[AX_APPEND_FLAG([--debug], [LDFLAGS_DEBUG])],
		[AC_MSG_WARN(['--debug' not supported. Are you sure this is an OmpSs compiler?])])
esac

AC_ARG_VAR([CFLAGS_PERF], [Extra CFLAGS for performance binary])
AC_ARG_VAR([LDFLAGS_PERF], [Extra LDFLAGS for performance binary])
AC_ARG_VAR([CFLAGS_INSTR], [Extra CFLAGS for instrumentation binary])
AC_ARG_VAR([LDFLAGS_INSTR], [Extra LDFLAGS for instrumentation binary])
AC_ARG_VAR([CFLAGS_DEBUG], [Extra CFLAGS for debug binary])
AC_ARG_VAR([LDFLAGS_DEBUG], [Extra LDFLAGS for debug binary])

AC_SUBST(CFLAGS_PERF)
AC_SUBST(LDFLAGS_PERF)
AC_SUBST(CFLAGS_INSTR)
AC_SUBST(LDFLAGS_INSTR)
AC_SUBST(CFLAGS_DEBUG)
AC_SUBST(LDFLAGS_DEBUG)

AC_ARG_ENABLE([sequence],
	AS_HELP_STRING([--disable-sequence], [disable sequence mode (overlapping of iterations); enabled by default]))

case "x${enable_ompss}x$enable_sequence" in
x*no*)
	AC_DEFINE([SEQUENCE_MODE], 0, [Disable sequence mode])
;;
*)
	AC_DEFINE([SEQUENCE_MODE], 1, [Enable sequence mode])
esac

AC_ARG_ENABLE([socket],
	AS_HELP_STRING([--disable-socket], [disable hard-coded socket-aware code; enabled by default for OmpSs]))

case "x${enable_ompss}x${enable_socket}" in
x*no*)
	AC_DEFINE([MANUAL_SOCKET], 0, [Disable hard-coded socket-aware])
;;
*)
	AC_DEFINE([MANUAL_SOCKET], 1, [Enable hard-coded socket-aware])
esac

AC_OUTPUT
