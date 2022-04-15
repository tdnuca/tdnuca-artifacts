MCC = smpcc
MLD = smpcc
LD = gcc
CC = gcc
CPPFLAGS = -I.
CFLAGS = -O3 -march=native -Wall -Werror -g
CFLAGS_seq = $(CFLAGS) -Wno-unknown-pragmas -Wno-unused-variable
CFLAGS_ompss = $(CFLAGS) --ompss --no-copy-deps
CFLAGS_instr = $(CFLAGS_ompss) --instrument
LDFLAGS = -g
LDFLAGS_seq = $(LDFLAGS)
LDFLAGS_ompss = $(LDFLAGS) --ompss --no-copy-deps
LDFLAGS_instr = $(LDFLAGS_ompss) --instrument
LDLIBS = -lm -lrt
LDLIBS_seq = $(LDLIBS)
LDLIBS_ompss = $(LDLIBS)
LDLIBS_instr = $(LDLIBS_ompss)

OBJS_seq = block-seq.o check_sum-seq.o comm-seq.o comm_alt-seq.o	\
	driver-seq.o init-seq.o main-seq.o move-seq.o plot-seq.o	\
	profile-seq.o refine-seq.o stencil-seq.o target-seq.o		\
	util-seq.o task_functions-seq.o

OBJS_ompss = block.o check_sum.o comm.o comm_alt.o driver.o init.o	\
	main.o move-seq.o plot-seq.o profile-seq.o refine-seq.o		\
	stencil.o target-seq.o util-seq.o task_functions-seq.o

OBJS_instr = block-i.o check_sum-i.o comm-i.o comm_alt-i.o driver-i.o	\
	init-i.o main-i.o move-seq.o plot-seq.o profile-seq.o		\
	refine-seq.o stencil-i.o target-seq.o util-seq.o		\
	task_functions-seq.o

HEADERS = block.h proto.h proto_task.h timer.h

all: ompss instr seq

seq: miniAMR-seq.x
ompss: miniAMR.x
instr: miniAMR-i.x

miniAMR-seq.x: $(OBJS_seq)
	$(LD) $(LDFLAGS_seq) -o $@ $^ $(LDLIBS_seq)

miniAMR.x: $(OBJS_ompss)
	$(MLD) $(LDFLAGS_ompss) -o $@ $^ $(LDLIBS_ompss)

miniAMR-i.x: $(OBJS_instr)
	$(MLD) $(LDFLAGS_instr) -o $@ $^ $(LDLIBS_instr)

%-seq.o: %.c $(HEADERS)
	$(CC) $(CPPFLAGS) $(CFLAGS_seq) -c -o $@ $<

%-i.o: %.c $(HEADERS)
	$(MCC) $(CPPFLAGS) $(CFLAGS_instr) -c -o $@ $<

%.o: %.c $(HEADERS)
	$(MCC) $(CPPFLAGS) $(CFLAGS_ompss) -c -o $@ $<

clean:
	rm -v -f *.o *~

distclean: clean
	rm -v -f miniAMR.x miniAMR-seq.x miniAMR-i.x

.PHONY: clean distclean seq ompss instr all


# Local variables:
# mode: makefile-gmake
# End:
