
CC=gcc
MCC=mcc
CFLAGS=-Wall  -O3 -DNDEBUG -fno-optimize-sibling-calls -fno-signaling-nans -ffast-math -fopenmp -Wno-error -Wno-unused-result -Wno-unused-variable -Wno-unused-local-typedefs -Wno-parentheses
OUTPUT=specfem3D

all:
	$(MCC) --ompss $(CFLAGS) specfem3D_smpss_v0.81.c  -o $(OUTPUT) -lpthread -lm

instr:
	$(MCC) --ompss $(CFLAGS) -instrument specfem3D_smpss_v0.81.c -o $(OUTPUT)_instr -lpthread -lm

clean:
	rm -f $(OUTPUT) $(OUTPUT)_instr *.o

.PHONY: all clean instr
