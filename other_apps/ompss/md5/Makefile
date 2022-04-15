
CC=gcc
MCC=mcc
CFLAGS=-Wall  -O3 -DNDEBUG -fno-optimize-sibling-calls -Wno-unused-variable -Wno-unused-local-typedefs -Wno-parentheses
LDFLAGS=-lpthread -lm
OBJS=md5_bmark.o md5.o
OUTPUT=md5

%.o: %.c
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
