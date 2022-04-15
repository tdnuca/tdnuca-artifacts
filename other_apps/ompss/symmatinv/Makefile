
CC=gcc
MCC=mcc
CFLAGS=--no-copy-deps -O3 -DDOUBLE_PREC -k -debug -g2 -frecord-gcc-switches --I${MKLROOT}/include/ -L${MKLROOT}/lib/intel64 -lmkl_sequential -lmkl_core -lmkl_rt
OBJS=inverse.o
OUTPUT=inverse

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
