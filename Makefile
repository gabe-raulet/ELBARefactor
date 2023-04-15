K?=7
MPICH=/usr/local/Cellar/mpich/4.1.1
MPICH_INC=-I$(MPICH)/include
MPICH_LIB=-L$(MPICH)/lib
COMPILER=g++-12
MPICH_FLAGS=$(MPICH_LIB) -L/usr/local/opt/libevent/lib -lmpi
FLAGS=-O0 -Wno-maybe-uninitialized -std=c++17 -DKMER_SIZE=$(K) $(MPICH_INC)

all: main

main: main.o CommGrid.o HashFuncs.o FastaIndex.o KmerComm.o Bloom.o Buffer.o HyperLogLog.o
	$(COMPILER) -o $@ $^ $(FLAGS) $(MPICH_FLAGS) -lz

test: test.o
	$(COMPILER) -o $@ $^ $(FLAGS) $(MPICH_FLAGS) -lz

main.o: main.cpp CommGrid.h Kmer.cpp Kmer.h
	$(COMPILER) $(FLAGS) -c -o $@ $<

test.o: test.cpp
	$(COMPILER) $(FLAGS) -c -o $@ $<

CommGrid.o: CommGrid.cpp CommGrid.h
	$(COMPILER) $(FLAGS) -c -o $@ $<

FastaIndex.o: FastaIndex.cpp FastaIndex.h CommGrid.h
	$(COMPILER) $(FLAGS) -c -o $@ $<

HyperLogLog.o: HyperLogLog.cpp HyperLogLog.h
	$(COMPILER) $(FLAGS) -c -o $@ $<

HashFuncs.o: HashFuncs.cpp HashFuncs.h
	$(COMPILER) $(FLAGS) -c -o $@ $<

KmerOps.o: KmerOps.cpp KmerOps.h
	$(COMPILER) $(FLAGS) -c -o $@ $<

KmerComm.o: KmerComm.cpp KmerComm.h
	$(COMPILER) $(FLAGS) -c -o $@ $<

Bloom.o: Bloom.cpp Bloom.h
	$(COMPILER) $(FLAGS) -c -o $@ $<

Buffer.o: Buffer.c Buffer.h
	gcc-12 -O2 -c -o $@ $<

clean:
	rm -rf *.o *.dSYM *.out
