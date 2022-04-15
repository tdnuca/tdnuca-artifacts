# target binaries
bin_PROGRAMS=cg cg_instr cg_seq


# A word on mcc and this file's structure

# For ompss targets, all non-ompss files will be compiled together in lib${target}.a with the backend compiler,
# while the other files will be compiled with mcc normally.
# This is so they can be given different compiler flags (without listing each .c -> .o individually).
# Hence target dependencies: binaries -> {ompss-only sources, lib.a -> {non-ompss sources}}

# At link-time mcc behaves like its backend with extra flags (-lnanos*, -L/path/to/nanos), no extra processing.


#define compiler/linker flags, common to all and specifics for ompss or plain C files
all_CFLAGS=-std=gnu11 -Wall -Wextra -Wpedantic -Werror -Wno-unused-result
all_CPPFLAGS=
all_LDFLAGS=

# add no-error for the ugliness mcc generates
ompss_CPPFLAGS=$(all_CPPFLAGS) -D_OMPSS
plain_CPPFLAGS=$(all_CPPFLAGS)

ompss_CFLAGS=$(all_CFLAGS) --ompss --Wn,-Wno-error=unused-result --Wn,-Wno-error=unused-parameter --Wn,-Wno-error=unused-but-set-variable --Wn,-Wno-error=unused-variable --Wn,-Wno-error=strict-aliasing --Wn,-Wno-error=pedantic
plain_CFLAGS=$(all_CFLAGS) -Wno-unknown-pragmas

ompss_LDFLAGS=$(all_LDFLAGS)
plain_LDFLAGS=$(all_LDFLAGS)

# explicitely pick compiler variant/flag when mcc is overkill
libcg_a-% libcg_instr_a-% cg_seq:CC=$(CC_BACKEND)
libcg_a-% libcg_instr_a-% cg_seq:CFLAGS=$(CFLAGS_BACKEND) $(all_CFLAGS)
libcg_a-% libcg_instr_a-% cg_seq:CPPFLAGS=$(CPPFLAGS_BACKEND) $(all_CPPFLAGS)
libcg_a-% libcg_instr_a-% cg_seq:LDFLAGS=$(LDFLAGS_BACKEND) $(all_LDFLAGS)

# target 'cg'
libcg_a_SOURCES=main.c matrix.c mmio.c debug.h global.h matrix.h mmio.h
libcg_a_CPPFLAGS=$(CPPFLAGS_BACKEND) $(plain_CPPFLAGS) -D_OMPSS
cg_SOURCES=cg.c counters.c cg.h counters.h
cg_CPPFLAGS=$(ompss_CPPFLAGS)
cg_CFLAGS=$(ompss_CFLAGS)
cg_LDFLAGS=$(ompss_LDFLAGS)
cg_LDADD=libcg.a


# target 'cg_instr' same as above with instrumentation flags
libcg_instr_a_SOURCES=main.c matrix.c mmio.c debug.h global.h matrix.h mmio.h
libcg_instr_a_CPPFLAGS=$(CPPFLAGS_BACKEND) $(plain_CPPFLAGS) -D_OMPSS
cg_instr_SOURCES=cg.c counters.c cg.h counters.h
cg_instr_CPPFLAGS=$(ompss_CPPFLAGS)
cg_instr_CFLAGS=$(ompss_CFLAGS) --instrument
cg_instr_LDFLAGS=$(ompss_LDFLAGS) --instrument
cg_instr_LDADD=libcg_instr.a


# target 'cg_seq', sequential baseline where all files are plain C and ompss pragmas are ignored
cg_seq_SOURCES=cg.c counters.c main.c matrix.c mmio.c cg.h counters.h debug.h global.h matrix.h mmio.h
cg_seq_CPPFLAGS=$(plain_CPPFLAGS)
cg_seq_CFLAGS=$(plain_CFLAGS)
cg_seq_LDFLAGS=$(plain_LDFLAGS)


# Do not install the lib*.a helpers
noinst_LIBRARIES=libcg.a libcg_instr.a


# if generating custom Extrae events is enabled, add necessary compiler & linkage flags
if USE_EXTRAE
cg_LDADD      +=$(LIBEXTRAE_NANOS)
cg_instr_LDADD+=$(LIBEXTRAE_NANOS)
cg_seq_LDADD   =$(LIBEXTRAE_SEQ)
endif

