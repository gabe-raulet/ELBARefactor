K?=31
L?=20
U?=30
BF?=1
COMPILE_TIME_PARAMETERS=-DKMER_SIZE=$(K) -DLOWER_KMER_FREQ=$(L) -DUPPER_KMER_FREQ=$(U) -DUSE_BLOOM=$(BF)
MPICH=/usr/local/Cellar/mpich/4.1.1
MPICH_INC=-I$(MPICH)/include
MPICH_LIB=-L$(MPICH)/lib
MPICH_FLAGS=
FLAGS=$(COMPILE_TIME_PARAMETERS) -O2 -Wno-maybe-uninitialized -Wno-deprecated -std=c++17 -I./inc -I./src

COMBBLAS=./CombBLAS
COMBBLAS_INC=$(COMBBLAS)/include/CombBLAS
COMBBLAS_SRC=$(COMBBLAS)/src
INCADD=-I$(COMBBLAS)/include/ -I$(COMBBLAS)/psort-1.0/include/ -I$(COMBBLAS)/usort/include/ -I$(COMBBLAS)/graph500-1.2/generator/include/

UNAME_S:=$(shell uname -s)

ifeq ($(UNAME_S),Linux)
COMPILER=mpic++
else ifeq ($(UNAME_S),Darwin)
COMPILER=g++-12
FLAGS+=$(MPICH_INC)
MPICH_FLAGS+=$(MPICH_LIB) -L/usr/local/opt/libevent/lib -lmpi
endif

all: elba

elba: main.o HashFuncs.o FastaIndex.o KmerComm.o Bloom.o HyperLogLog.o ReadOverlap.o CommGrid.o MPIType.o Logger.o
	@echo CXX -c -o $@ $^
	@$(COMPILER) $(FLAGS) $(INCADD) -o $@ $^ $(MPICH_FLAGS) -lz

%.o: src/%.cpp
	@echo CXX $(COMPILE_TIME_PARAMETERS) -c -o $@ $<
	@$(COMPILER) $(FLAGS) $(INCADD) -c -o $@ $<

main.o: src/main.cpp src/Kmer.cpp inc/Kmer.h
FastaIndex.o: src/FastaIndex.cpp inc/FastaIndex.h
HyperLogLog.o: src/HyperLogLog.cpp inc/HyperLogLog.h
HashFuncs.o: src/HashFuncs.cpp inc/HashFuncs.h
KmerComm.o: src/KmerComm.cpp inc/KmerComm.h inc/Bloom.h
Bloom.o: src/Bloom.cpp inc/Bloom.h
ReadOverlap.o: src/ReadOverlap.cpp inc/ReadOverlap.h
Logger.o: src/Logger.cpp inc/Logger.h

CommGrid.o: $(COMBBLAS_SRC)/CommGrid.cpp $(COMBBLAS_INC)/CommGrid.h
	@echo CXX -c -o $@ $<
	@$(COMPILER) $(FLAGS) $(INCADD) -c -o $@ $<

MPIType.o: $(COMBBLAS_SRC)/MPIType.cpp $(COMBBLAS_INC)/MPIType.h
	@echo CXX -c -o $@ $<
	@$(COMPILER) $(FLAGS) $(INCADD) -c -o $@ $<

clean:
	rm -rf *.o *.dSYM *.out

gitclean: clean
	git clean -f
