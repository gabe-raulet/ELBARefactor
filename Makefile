MPICH=/usr/local/Cellar/mpich/4.1.1
MPICH_INC=-I$(MPICH)/include
MPICH_LIB=-L$(MPICH)/lib
COMPILER=g++-12
ALL_FLAGS=-O2 -Wno-maybe-uninitialized -std=c++17 $(MPICH_INC) $(MPICH_LIB) -L/usr/local/opt/libevent/lib -lmpi
FLAGS=-O2 -Wno-maybe-uninitialized -std=c++17 $(MPICH_INC)

all: main

main: main.o CommGrid.o HashFuncs.o FastaIndex.o KmerOps.o Bloom.o Buffer.o HyperLogLog.o
	$(COMPILER) -o $@ $^ $(ALL_FLAGS) -lz

main.o: main.cpp CommGrid.h Kmer.cpp Kmer.h
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

Bloom.o: Bloom.cpp Bloom.h
	$(COMPILER) $(FLAGS) -c -o $@ $<

Buffer.o: Buffer.c Buffer.h
	gcc-12 -O2 -c -o $@ $<

clean:
	rm -rf *.o *.dSYM *.out
