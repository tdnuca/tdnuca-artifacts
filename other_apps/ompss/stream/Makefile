CC = smpcc
OPT=3
CFLAGS  = -O$(OPT) -ffast-math $(OMPSSFLAGS) -Wall -Wextra
OMPSSFLAGS = --ompss
CPPFLAGS =
LDFLAGS =
LDLIBS = -lm
SOURCES = stream.c
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = stream
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
	$(CC) $(CFLAGS) $(LDFLAGS) $< -o $@ $(LDLIBS)

%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $< -c -o $@

clean:
	rm -rf *.o $(EXECUTABLES)

.INTERMEDIATE: $(SOURCES:.c=.o)
.PHONY: perf instr
