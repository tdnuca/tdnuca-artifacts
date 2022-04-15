# SYNOPSIS
#
#   forked from AX_PROG_CC_FOR_BUILD
#   http://www.gnu.org/software/autoconf-archive/ax_prog_cc_for_build.html
#
# DESCRIPTION
#
#   This macro searches for a C compiler that generates native executables,
#   that is a C compiler that surely is not a cross-compiler. This can be
#   useful if you have to generate source code at compile-time like for
#   example GCC does.
#
#   The macro sets the CC_FOR_OMPSS and CPP_FOR_OMPSS macros to anything
#   needed to compile or link (CC_FOR_OMPSS) and preprocess (CPP_FOR_OMPSS).
#   The value of these variables can be overridden by the user by specifying
#   a compiler with an environment variable (like you do for standard CC).
#
#   It also sets BUILD_EXEEXT and BUILD_OBJEXT to the executable and object
#   file extensions for the build platform, and GCC_FOR_OMPSS to `yes' if
#   the compiler we found is GCC. All these variables but GCC_FOR_OMPSS are
#   substituted in the Makefile.
#
# LICENSE
#
#   Copyright (c) 2008 Paolo Bonzini <bonzini@gnu.org>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 8


AU_ALIAS([AC_PROG_CC_FOR_OMPSS], [AX_PROG_CC_FOR_OMPSS])
AC_DEFUN([AX_PROG_CC_FOR_OMPSS], [dnl
AC_REQUIRE([AC_PROG_CC])dnl
AC_REQUIRE([AC_PROG_CPP])dnl
AC_REQUIRE([AC_EXEEXT])dnl
AC_REQUIRE([AC_CANONICAL_HOST])dnl

m4_version_prereq([2.69],[
	dnl hack to workaround the fact that push/popdef doesn't propagate anymore
	dnl inside functions like ac_fn_c_try_link (ie expanding old values, after push)
	CC_SAVED="$CC"
	CPP_SAVED="$CPP"
	CFLAGS_SAVED="$CFLAGS"
	CPPFLAGS_SAVED="$CPPFLAGS"
	LDFLAGS_SAVED="$LDFLAGS"

	export CC="${CC_FOR_OMPSS:-$CC}"
	export CPP="${CPP_FOR_OMPSS:-$CPP}"
	export CFLAGS="${CFLAGS_FOR_OMPSS:-$CFLAGS}"
	export CPPFLAGS="${CPPFLAGS_FOR_OMPSS:-$CPPFLAGS}"
	export LDFLAGS="${LDFLAGS_FOR_OMPSS:-$LDFLAGS}"
],[
])

dnl variants of CPPFLAGS
export CPPFLAGS_FOR_OMPSS_BARE="$CPPFLAGS_FOR_OMPSS"

dnl Use the standard macros, but make them use other variable names
dnl
pushdef([ac_cv_prog_CPP], ac_cv_ompss_prog_CPP)dnl
pushdef([ac_cv_prog_gcc], ac_cv_ompss_prog_gcc)dnl
pushdef([ac_cv_prog_cc_works], ac_cv_ompss_prog_cc_works)dnl
pushdef([ac_cv_prog_cc_cross], ac_cv_ompss_prog_cc_cross)dnl
pushdef([ac_cv_prog_cc_g], ac_cv_ompss_prog_cc_g)dnl
pushdef([ac_cv_exeext], ac_cv_ompss_exeext)dnl
pushdef([ac_cv_objext], ac_cv_ompss_objext)dnl
pushdef([ac_exeext], ac_ompss_exeext)dnl
pushdef([ac_objext], ac_ompss_objext)dnl
pushdef([CC], CC_FOR_OMPSS)dnl
pushdef([CPP], CPP_FOR_OMPSS)dnl
pushdef([CFLAGS], CFLAGS_FOR_OMPSS)dnl
pushdef([LDFLAGS], LDFLAGS_FOR_OMPSS)dnl
pushdef([CPPFLAGS], CPPFLAGS_FOR_OMPSS)dnl
pushdef([host], build)dnl
pushdef([host_alias], build_alias)dnl
pushdef([host_cpu], build_cpu)dnl
pushdef([host_vendor], build_vendor)dnl
pushdef([host_os], build_os)dnl
pushdef([ac_cv_host], ac_cv_ompss)dnl
pushdef([ac_cv_host_alias], ac_cv_ompss_alias)dnl
pushdef([ac_cv_host_cpu], ac_cv_ompss_cpu)dnl
pushdef([ac_cv_host_vendor], ac_cv_ompss_vendor)dnl
pushdef([ac_cv_host_os], ac_cv_ompss_os)dnl
pushdef([ac_cpp], ac_ompss_cpp)dnl
pushdef([ac_compile], ac_ompss_compile)dnl
pushdef([ac_link], ac_ompss_link)dnl


save_cross_compiling=$cross_compiling
save_ac_tool_prefix=$ac_tool_prefix
cross_compiling=no
ac_tool_prefix=

AC_PROG_CC
AC_PROG_CPP
AC_EXEEXT


m4_version_prereq([2.69],[
	dnl hack to workaround the fact that push/popdef doesn't propagate anymore
	dnl inside functions like ac_fn_c_try_link (ie expanding old values, after push)
	popdef([CPPFLAGS])
	export CPPFLAGS="$CPPFLAGS_FOR_OMPSS_BARE --do-not-process-file"
	pushdef([CPPFLAGS], CPPFLAGS_FOR_OMPSS)
],[
])
export CPPFLAGS_FOR_OMPSS="$CPPFLAGS_FOR_OMPSS_BARE --do-not-process-file"

AC_MSG_CHECKING([whether $CC can be bypassed to compile C code (using --do-not-process-file)])
AC_LINK_IFELSE(
	[AC_LANG_PROGRAM()],
	[AC_MSG_RESULT([yes])], [AC_MSG_FAILURE([$CC can't compile normal C code])]
)

m4_version_prereq([2.69],[
	dnl hack to workaround the fact that push/popdef doesn't propagate anymore
	dnl inside functions like ac_fn_c_try_link (ie expanding old values, after push)
	popdef([CPPFLAGS])
	export CPPFLAGS="$CPPFLAGS_FOR_OMPSS_BARE --ompss"
	pushdef([CPPFLAGS], CPPFLAGS_FOR_OMPSS)
],[
])
export CPPFLAGS_FOR_OMPSS="$CPPFLAGS_FOR_OMPSS_BARE --ompss"

AC_MSG_CHECKING([whether $CC can compile OmpSs code (using --ompss)])
AC_LINK_IFELSE(
	[AC_LANG_PROGRAM(
		[[
		#include <stdio.h>
		#include <nanox/nanos_omp.h>
		]],
		[[
		char str[] = "Hello World on %d threads!\n";
		// making it clearly OmpSs not OpenMP dependencies
		#pragma omp task inout(str)
		{
			printf(str, nanos_omp_get_num_threads());
		}
		#pragma omp taskwait on(str)
		]]
	)],
	[AC_MSG_RESULT([yes])], [AC_MSG_FAILURE([$CC can't compile C code with OmpSs pragmas])]
)

m4_version_prereq([2.69],[
	dnl hack to workaround the fact that push/popdef doesn't propagate anymore
	dnl inside functions like ac_fn_c_try_link (ie expanding old values, after push)
	popdef([CPPFLAGS])
	export CPPFLAGS="$CPPFLAGS_FOR_OMPSS_BARE"
	export CPPFLAGS_FOR_OMPSS="$CPPFLAGS_FOR_OMPSS_BARE"
	pushdef([CPPFLAGS], CPPFLAGS_FOR_OMPSS)
],[
])
export CPPFLAGS_FOR_OMPSS="$CPPFLAGS_FOR_OMPSS_BARE"

ac_tool_prefix=$save_ac_tool_prefix
cross_compiling=$save_cross_compiling

dnl Restore the old definitions

dnl
popdef([ac_link])dnl
popdef([ac_compile])dnl
popdef([ac_cpp])dnl
popdef([ac_cv_host_os])dnl
popdef([ac_cv_host_vendor])dnl
popdef([ac_cv_host_cpu])dnl
popdef([ac_cv_host_alias])dnl
popdef([ac_cv_host])dnl
popdef([host_os])dnl
popdef([host_vendor])dnl
popdef([host_cpu])dnl
popdef([host_alias])dnl
popdef([host])dnl
popdef([CPPFLAGS])
popdef([LDFLAGS])dnl
popdef([CFLAGS])dnl
popdef([CPP])dnl
popdef([CC])dnl
popdef([ac_objext])dnl
popdef([ac_exeext])dnl
popdef([ac_cv_objext])dnl
popdef([ac_cv_exeext])dnl
popdef([ac_cv_prog_cc_g])dnl
popdef([ac_cv_prog_cc_cross])dnl
popdef([ac_cv_prog_cc_works])dnl
popdef([ac_cv_prog_gcc])dnl
popdef([ac_cv_prog_CPP])dnl

dnl Finally, set Makefile variables
dnl
BUILD_EXEEXT=$ac_ompss_exeext
BUILD_OBJEXT=$ac_ompss_objext
AC_SUBST(BUILD_EXEEXT)dnl
AC_SUBST(BUILD_OBJEXT)dnl
AC_SUBST([CFLAGS_FOR_OMPSS])dnl
AC_SUBST([CPPFLAGS_FOR_OMPSS])dnl
AC_SUBST([LDFLAGS_FOR_OMPSS])dnl

m4_version_prereq([2.69],[
	dnl hack to workaround the fact that push/popdef doesn't propagate anymore
	dnl inside functions like ac_fn_c_try_link (ie expanding old values, after push)
	export CC="$CC_SAVED"
	export CPP="$CPP_SAVED"
	export CFLAGS="$CFLAGS_SAVED"
	export CPPFLAGS="$CPPFLAGS_SAVED"
	export LDFLAGS="$LDFLAGS_SAVED"
],[
])

])

