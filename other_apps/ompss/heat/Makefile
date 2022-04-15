CXX=g++
MCXX=smpcxx
CXXFLAGS+=-Wall -std=c++03 -O3 $(LDFLAGS)
CXXFLAGS_seq=-Wno-unknown-pragmas
CXXFLAGS_debug=-O0 -Wall -std=c++03 -g -ggdb3
CXXFLAGS_ompss=--ompss -Wno-unused-local-typedefs

%-seq.o: %.cc
	$(CXX) $(CXXFLAGS) $(CXXFLAGS_seq) -c -o $@ $^

%.o: %.cc
	$(MCXX) $(CXXFLAGS) $(CXXFLAGS_ompss) --no-copy-deps -c -o $@ $^

%-cd.o: %.cc
	$(MCXX) $(CXXFLAGS) $(CXXFLAGS_ompss) --copy-deps -c -o $@ $^

%-i.o: %.cc
	$(MCXX) $(CXXFLAGS) $(CXXFLAGS_ompss) --no-copy-deps --instrument -c -o $@ $^

%-d.o: %.cc
	$(MCXX) $(CXXFLAGS_debug) $(CXXFLAGS_ompss) --no-copy-deps --debug -k -c -o $@ $^

%-seq-d.o: %.cc
	$(CXX) $(CXXFLAGS_debug) $(CXXFLAGS_seq) -c -o $@ $^

%-i-cd.o: %.cc
	$(MCXX) $(CXXFLAGS) $(CXXFLAGS_ompss) --copy-deps --instrument -c -o $@ $^

%-d-cd.o: %.cc
	$(MCXX) $(CXXFLAGS_debug) $(CXXFLAGS_ompss) --copy-deps --debug -k -c -o $@ $^

%-noinit.o: %.cc
	$(MCXX) $(CXXFLAGS) $(CXXFLAGS_ompss) --no-copy-deps -DINIT_TASKS=0 -c -o $@ $^

%-i-noinit.o: %.cc
	$(MCXX) $(CXXFLAGS) $(CXXFLAGS_ompss) --no-copy-deps --instrumentation -DINIT_TASKS=0 -c -o $@ $^

%-d-noinit.o: %.cc
	$(MCXX) $(CXXFLAGS_debug) $(CXXFLAGS_ompss) --no-copy-deps --debug -DINIT_TASKS=0 -c -o $@ $^

%-noinit-cd.o: %.cc
	$(MCXX) $(CXXFLAGS) $(CXXFLAGS_ompss) --copy-deps -DINIT_TASKS=0 -c -o $@ $^

%-i-noinit-cd.o: %.cc
	$(MCXX) $(CXXFLAGS) $(CXXFLAGS_ompss) --copy-deps --instrumentation -DINIT_TASKS=0 -c -o $@ $^

%-d-noinit-cd.o: %.cc
	$(MCXX) $(CXXFLAGS_ompss) $(CXXFLAGS_debug) --copy-deps --debug -DINIT_TASKS=0 -c -o $@ $^

heat-seq: heat-seq.o matrix-seq.o misc-seq.o algorithms-seq.o algorithms_inner-seq.o
	$(CXX) $(CXXFLAGS) $(CXXFLAGS_seq) -o $@ $^ $(LIBS)

heat: heat.o matrix.o misc-seq.o algorithms.o algorithms_inner-seq.o
	$(MCXX) $(CXXFLAGS) $(CXXFLAGS_ompss) -o $@ $^ $(LIBS)

heat-noinit: heat.o matrix-noinit.o misc-seq.o algorithms.o algorithms_inner-seq.o
	$(MCXX) $(CXXFLAGS) $(CXXFLAGS_ompss) -DINIT_TASKS=0 -o $@ $^ $(LIBS)

heat-i-noinit: heat-i.o matrix-i-noinit.o misc-seq.o algorithms-i.o algorithms_inner-seq.o
	$(MCXX) $(CXXFLAGS) $(CXXFLAGS_ompss) --instrumentation -DINIT_TASKS=0 -o $@ $^ $(LIBS)

heat-i: heat-i.o matrix-i.o misc-seq.o algorithms-i.o algorithms_inner-seq.o
	$(MCXX) $(CXXFLAGS) $(CXXFLAGS_ompss) --instrument -o $@ $^ $(LIBS)

heat-d: heat-d.o matrix-d.o misc-seq-d.o algorithms-d.o algorithms_inner-seq-d.o
	$(MCXX) $(CXXFLAGS_debug) $(CXXFLAGS_ompss) --debug -o $@ $^ $(LIBS)

heat-seq-d: heat-seq-d.o matrix-seq-d.o misc-seq-d.o algorithms-seq-d.o algorithms_inner-seq-d.o
	$(CXX) $(CXXFLAGS_debug) $(CXXFLAGS_seq) -o $@ $^ $(LIBS)

heat-cd: heat-cd.o matrix-cd.o misc-seq.o algorithms-cd.o algorithms_inner-seq.o
	$(MCXX) $(CXXFLAGS) $(CXXFLAGS_ompss) -o $@ $^ $(LIBS)

heat-noinit-cd: heat-cd.o matrix-noinit-cd.o misc-seq.o algorithms-cd.o algorithms_inner-seq.o
	$(MCXX) $(CXXFLAGS) $(CXXFLAGS_ompss) -DINIT_TASKS=0 -o $@ $^ $(LIBS)

heat-i-noinit-cd: heat-i.o matrix-i-noinit-cd.o misc-seq.o algorithms-i-cd.o algorithms_inner-seq.o
	$(MCXX) $(CXXFLAGS) $(CXXFLAGS_ompss) --instrumentation -DINIT_TASKS=0 -o $@ $^ $(LIBS)

heat-i-cd: heat-i-cd.o matrix-i-cd.o misc-seq.o algorithms-i-cd.o algorithms_inner-seq.o
	$(MCXX) $(CXXFLAGS) $(CXXFLAGS_ompss) --instrument -o $@ $^ $(LIBS)

heat-d-cd: heat-d-cd.o matrix-d-cd.o misc-seq-d.o algorithms-d-cd.o algorithms_inner-seq-d.o
	$(MCXX) $(CXXFLAGS_debug) $(CXXFLAGS_ompss) --debug -o $@ $^ $(LIBS)

all: heat heat-seq heat-i heat-d heat-seq-d heat-noinit heat-i-noinit heat-cd heat-i-cd heat-d-cd heat-noinit-cd heat-i-noinit-cd 

default: heat heat-cd

.PHONY: clean distclean all

.DEFAULT_GOAL=default

clean:
	rm -fv *.o *~

distclean: clean
	rm -fv heat heat-seq heat-i heat-d heat-seq-d heat-noinit heat-i-noinit heat-cd heat-seq-cd heat-seq-d-cd heat-i-cd heat-d-cd heat-noinit-cd heat-i-noinit-cd $(MCXX)_*.cc

install:
	mkdir -p bin
	mv heat bin/heat-ompss || :
	mv heat-seq bin/heat-serial || :
