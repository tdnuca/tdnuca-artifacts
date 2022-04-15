MCC=smpcc
MCXX=smpcxx
CXX=g++
FC=gfortran

CXXFLAGS=-Wall -O3 -march=native -std=c++03
CXXFLAGS_seq=-Wno-unknown-pragmas
CXXFLAGS_debug=-g -ggdb3
CXXFLAGS_ompss=--ompss -Wno-unused-local-typedefs

LIBS=-L$(MY_SCRATCH)/compiled/openblas/lib -lopenblas -lgfortran

%-seq.o: %.cc Makefile
	$(CXX) $(CXXFLAGS) $(CXXFLAGS_seq) -c -o $@ $<

blocksize.o: blocksize.f90 Makefile
	$(FC) -c -o $@ $<

%.o: %.cc Makefile
	$(MCXX) $(CXXFLAGS) $(CXXFLAGS_ompss) -c -o $@ $<

%-i.o: %.cc Makefile
	$(MCXX) $(CXXFLAGS) $(CXXFLAGS_ompss) --instrument -c -o $@ $<

%-d.o: %.cc Makefile
	$(MCXX) $(CXXFLAGS) $(CXXFLAGS_debug) $(CXXFLAGS_ompss) --debug -k -c -o $@ $<

qr-seq: qr-tile-seq.o qr-seq.o misc-seq.o matrix-seq.o blocksize.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

qr: qr-tile.o qr.o misc-seq.o matrix.o blocksize.o
	$(MCXX) $(CXXFLAGS) $(CXXFLAGS_ompss) -o $@ $^ $(LIBS)

qr-i: qr-tile-i.o qr-i.o misc-seq.o matrix-i.o blocksize.o
	$(MCXX) $(CXXFLAGS) $(CXXFLAGS_ompss) --instrument -o $@ $^ $(LIBS)

qr-d: qr-tile-d.o qr-d.o misc-seq.o matrix-d.o blocksize.o
	$(MCXX) $(CXXFLAGS) $(CXXFLAGS_ompss) $(CXXFLAGS_debug) --debug -o $@ $^ $(LIBS)

all: qr qr-seq qr-i qr-d

.PHONY: clean distclean all

.DEFAULT_GOAL=all

clean:
	rm -fv *.o *~

distclean: clean
	rm -fv qr qr-seq qr-i qr-d $(MCXX)_*.cc
