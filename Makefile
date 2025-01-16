# for linux:
# apt install armadillo
# for macos:
# brew install armadillo
CXX=g++
CC=g++ # use the C++ compiler for linking

# the default Apple clang doesn't support OpenMP
ifeq ($(shell uname), Darwin)
CXX=g++-14
CC=g++-14
endif

COMP_FLAGS=-O3 -Wall -std=c++14
ARMA_FLAGS=$(shell pkg-config --cflags armadillo)
ARMA_LIB=$(shell pkg-config --libs armadillo)
OMP_FLAGS=-fopenmp

all: build solve simpsolve

build: CXXFLAGS = $(COMP_FLAGS) $(OMP_FLAGS)
build: CPPFLAGS = $(ARMA_FLAGS)
build: LDFLAGS = $(OMP_FLAGS) $(ARMA_LIB)
build: build.o RWGDomain.o fio.o

solve: CXXFLAGS = $(COMP_FLAGS)
solve: CPPFLAGS = $(ARMA_FLAGS)
solve: LDLIBS = $(ARMA_LIB)
solve: solve.o RWGDomain.o fio.o

simpsolve: CXXFLAGS = $(COMP_FLAGS)
simpsolve: CPPFLAGS = $(ARMA_FLAGS)
simpsolve: LDLIBS = $(ARMA_LIB)
simpsolve: simpsolve.o fio.o

integrate.o: point.h integrate.h
RWGDomain.o: RWGDomain.cc point.h RWGDomain.h

silo: data2_current.silo data3_current.silo data4_current.silo data5_current.silo
data2_current.silo: data2.tri data2.curJ
	-emsurftranslator -s data2
data3_current.silo: data3.tri data3.curJ
	-emsurftranslator -s data3
data4_current.silo: data4.tri data4.curJ
	-emsurftranslator -s data4
data5_current.silo: data5.tri data5.curJ
	-emsurftranslator -s data5

clean:
	rm -f *.o proj *.gnu build solve simpsolve

