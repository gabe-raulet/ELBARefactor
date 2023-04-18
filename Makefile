K?=31
L?=20
U?=30
BF?=1
COMPILE_TIME_PARAMETERS=-DKMER_SIZE=$(K) -DLOWER_KMER_FREQ=$(L) -DUPPER_KMER_FREQ=$(U) -DUSE_BLOOM=$(BF)
MPICH=/usr/local/Cellar/mpich/4.1.1
MPICH_INC=-I$(MPICH)/include
MPICH_LIB=-L$(MPICH)/lib
MPICH_FLAGS=
FLAGS=$(COMPILE_TIME_PARAMETERS) -O2 -Wno-maybe-uninitialized -std=c++17
UNAME_S:=$(shell uname -s)

ifeq ($(UNAME_S),Linux)
COMPILER=mpic++
else ifeq ($(UNAME_S),Darwin)
COMPILER=g++-12
FLAGS+=$(MPICH_INC)
MPICH_FLAGS+=$(MPICH_LIB) -L/usr/local/opt/libevent/lib -lmpi
endif

all: main

main: main.o CommGrid.o HashFuncs.o FastaIndex.o KmerComm.o Bloom.o HyperLogLog.o
	$(COMPILER) -o $@ $^ $(FLAGS) $(MPICH_FLAGS) -lz

test: test.o
	$(COMPILER) -o $@ $^ $(FLAGS) $(MPICH_FLAGS) -lz

main.o: main.cpp CommGrid.h Kmer.cpp Kmer.h
	@echo CXX $(COMPILE_TIME_PARAMETERS) -c -o $@ $<
	@$(COMPILER) $(FLAGS) -c -o $@ $<

test.o: test.cpp
	@echo CXX $(COMPILE_TIME_PARAMETERS) -c -o $@ $<
	$(COMPILER) $(FLAGS) -c -o $@ $<

CommGrid.o: CommGrid.cpp CommGrid.h
	@echo CXX $(COMPILE_TIME_PARAMETERS) -c -o $@ $<
	@$(COMPILER) $(FLAGS) -c -o $@ $<

FastaIndex.o: FastaIndex.cpp FastaIndex.h CommGrid.h
	@echo CXX $(COMPILE_TIME_PARAMETERS) -c -o $@ $<
	@$(COMPILER) $(FLAGS) -c -o $@ $<

HyperLogLog.o: HyperLogLog.cpp HyperLogLog.h
	@echo CXX $(COMPILE_TIME_PARAMETERS) -c -o $@ $<
	@$(COMPILER) $(FLAGS) -c -o $@ $<

HashFuncs.o: HashFuncs.cpp HashFuncs.h
	@echo CXX $(COMPILE_TIME_PARAMETERS) -c -o $@ $<
	@$(COMPILER) $(FLAGS) -c -o $@ $<

KmerOps.o: KmerOps.cpp KmerOps.h
	@echo CXX $(COMPILE_TIME_PARAMETERS) -c -o $@ $<
	@$(COMPILER) $(FLAGS) -c -o $@ $<

KmerComm.o: KmerComm.cpp KmerComm.h Bloom.h
	@echo CXX $(COMPILE_TIME_PARAMETERS) -c -o $@ $<
	@$(COMPILER) $(FLAGS) -c -o $@ $<

Bloom.o: Bloom.cpp Bloom.h
	@echo CXX $(COMPILE_TIME_PARAMETERS) -c -o $@ $<
	@$(COMPILER) $(FLAGS) -c -o $@ $<

clean:
	rm -rf *.o *.dSYM *.out
