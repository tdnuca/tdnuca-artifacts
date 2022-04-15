EXEC_PRE=prk2_stencil

ifdef version
  ifeq "$(version)" "ompss"
    CC = mcc
    CXXFLAGS += -g -O3 --ompss -I. -Wall
    LDFLAGS += -g -O3 --ompss
    EXEC = $(EXEC_PRE)-ompss
  endif
  ifeq "$(version)" "serial"
    #Below are reasonable default flags for a serial build
    CC = mcc
    CXXFLAGS += -g -O3 --ompss -I. -Wall
    LDFLAGS += -g -O3 --ompss
    EXEC = $(EXEC_PRE)-serial
  endif
endif

all: $(EXEC)

clean:
	rm -f stencil_kernel.o prk2_stencil-*

stencil_kernel.o: src/stencil_kernel.c src/stencil_kernel.h
	$(CC) $(CXXFLAGS) -restrict -strict-ansi -c -o stencil_kernel.o src/stencil_kernel.c

$(EXEC): src/intel_prk2_stencil.c src/stencil_kernel.h stencil_kernel.o
	$(CC) $(CXXFLAGS) src/intel_prk2_stencil.c -o prk2_stencil-$(version) $(LDFLAGS) $(LIBS) stencil_kernel.o

install:
	mkdir -p ./bin
	mv $(EXEC_PRE)-* ./bin
