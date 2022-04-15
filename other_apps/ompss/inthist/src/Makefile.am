bin_PROGRAMS = inthist_crossweave inthist_crossweave_i inthist_crossweave_d
DISTCLEANFILES = *.in

ih_sources = inthist_check.c \
	inthist_kernels.c \
	inthist_randomimage.c \
	inthist_utils.c \
	inthist_main.c

ih_headers = inthist_extra.h \
	inthist_check.h \
	inthist_image.h \
	inthist_kernels.h \
	inthist_utils.h

inthist_crossweave_SOURCES = $(ih_sources) $(ih_headers)
inthist_crossweave_i_SOURCES = $(ih_sources) $(ih_headers)
inthist_crossweave_d_SOURCES = $(ih_sources) $(ih_headers)

inthist_crossweave_CFLAGS = $(CFLAGS_PERF)
inthist_crossweave_LDFLAGS = $(LDFLAGS_PERF)

inthist_crossweave_i_CFLAGS = $(CFLAGS_INSTR)
inthist_crossweave_i_LDFLAGS = $(LDFLAGS_INSTR)

inthist_crossweave_d_CFLAGS = $(CFLAGS_DEBUG)
inthist_crossweave_d_LDFLAGS = $(LDFLAGS_DEBUG)
