
CC=gcc
MCC=mcc
CFLAGS=-Wall  -O3 -DNDEBUG -Wno-unused-result -fno-optimize-sibling-calls -Wno-unused-variable -Wno-unused-local-typedefs -Wno-parentheses
OUTPUT=stap

all:
	$(MCC) --ompss $(CFLAGS) stap_bench_trans.c -o $(OUTPUT)

instr:
	$(MCC) --ompss $(CFLAGS) -instrument stap_bench_trans.c -o $(OUTPUT)_instr

clean:
	rm $(OUTPUT) $(OUTPUT)_instr

.PHONY: all clean instr
