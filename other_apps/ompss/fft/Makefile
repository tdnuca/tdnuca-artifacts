CC = smpcc
OPT=3
CFLAGS  = -O$(OPT) -ffast-math $(OMPSSFLAGS) -Wall -Wextra
OMPSSFLAGS = --ompss
CPPFLAGS = -D_Float128=__float128
LDFLAGS =
LIBS = -lm
SOURCES = fft.c seq.c benchmark.c
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = fft
EXECUTABLES = $(EXECUTABLE) $(EXECUTABLE)_instr $(EXECUTABLE)_debug $(EXECUTABLE)_seq


do: perf

perf: $(EXECUTABLE)

instr: OMPSSFLAGS += --instrument
instr: $(EXECUTABLE)_instr

debug: OMPSSFLAGS += --debug -K
debug: CFLAGS += -g
debug: OPT = 0
debug: $(EXECUTABLE)_debug

seq:$(EXECUTABLE)_seq
seq:CC = gcc
seq:OMPSSFLAGS =  -Wno-unknown-pragmas


$(EXECUTABLES): $(OBJECTS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(OBJECTS) -o $@ $(LIBS)

%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

clean:
	rm -rf *.o $(EXECUTABLES)

.INTERMEDIATE: $(SOURCES:.c=.o)
.PHONY: perf instr
