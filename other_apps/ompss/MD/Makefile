CC = gcc
MCC = smpcc --cc=$(CC)

CFLAGS=-O$(O) -std=gnu11 -rdynamic -Wall -Wextra -g -rdynamic
OMPSSFLAGS = --ompss -K -D_Float128=__float128 --no-copy-deps --Wn,-Wno-unused-but-set-variable,-Wno-unused-parameter
CPPFLAGS += -D_GNU_SOURCE
LDFLAGS=
LDLIBS=-lm

O=3
%_debug:O=0

# un-strict aliasing needed for the ways #atomic on doubles is translated into __sync_cas
%_perf.o: OMPSSFLAGS+=--keep-all-files --output-dir=.build_perf
%_debug.o:OMPSSFLAGS+=--keep-all-files --output-dir=.build_debug --debug
%_instr.o:OMPSSFLAGS+=--keep-all-files --output-dir=.build_instr --instrument

%_debug:OMPSSFLAGS+=--debug
%_instr:OMPSSFLAGS+=--instrument

md_seq md_seq.o: MCC=$(CC)
md_seq io_lammps.o: OMPSSFLAGS=
md_seq.o:OMPSSFLAGS=-Wno-unknown-pragmas -D_SEQUENTIAL

CFLAGS+=$(OMPSSFLAGS)

all:   perf instr debug seq
seq:   md_seq
perf:  md_perf
instr: md_instr
debug: md_debug

md_%.o: md.c | dirs
	$(MCC) $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -c $< -o $@

md_%: md_%.o io_lammps.o
	$(MCC) $(CFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

.PHONY: clean dirs all seq perf instr debug

dirs:
	@mkdir -p .build_perf .build_debug .build_instr

clean:
	@rm -rfv *.o md_seq md_perf md_instr md_debug
