
CC=g++
MCC=mcxx

CFLAGS=-Wall  -O3 -DNDEBUG -Wno-sign-compare -fno-optimize-sibling-calls -Wno-unused-variable -Wno-unused-local-typedefs -Wno-parentheses
LDFLAGS=-lpthread -lm
OBJS=benchmark_engine.o convert_engine.o image.o program.o rotation_engine.o
OUTPUT=rgbyuv

%.o: %.cpp
	$(MCC) --ompss $(CFLAGS) -c $< -o $@

all: $(OBJS)
	$(MCC) --ompss $(CFLAGS) $(OBJS) -o $(OUTPUT) $(LDFLAGS)
	rm -f $(OBJS)

instr: $(OBJS)
	$(MCC) --ompss $(CFLAGS) -instrument $(OBJS) -o $(OUTPUT)_instr $(LDFLAGS)
	rm -f $(OBJS)

clean:
	rm -f $(OUTPUT) $(OUTPUT)_instr *.o

.PHONY: all clean instr
