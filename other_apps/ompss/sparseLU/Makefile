
CC=gcc
MCC=mcc
CFLAGS=-Wall  -O3 -DNDEBUG -fno-optimize-sibling-calls -Wno-unused-variable -Wno-unused-local-typedefs -Wno-parentheses
OUTPUT=sparseLU

all:
	$(MCC) --ompss $(CFLAGS) checkSparseLU.c -o $(OUTPUT) -lpthread -lm

instr:
	$(MCC) --ompss $(CFLAGS) -instrument checkSparseLU.c -o $(OUTPUT)_instr -lpthread -lm

clean:
	rm -f $(OUTPUT) $(OUTPUT)_instr

.PHONY: all clean instr
