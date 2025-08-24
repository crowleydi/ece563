# for linux:
# apt install armadillo
# for macos:
# brew install armadillo
CXX=g++
CC=g++ # use the C++ compiler for linking

# the default Apple clang doesn't support OpenMP
ifeq ($(shell uname), Darwin)
CXX=g++-15
CC=g++-15
endif

COMP_FLAGS=-O3 -Wall -std=c++17
ARMA_FLAGS=$(shell pkg-config --cflags armadillo)
ARMA_LIB=$(shell pkg-config --libs armadillo)
OMP_FLAGS=-fopenmp

all: build solve

build: CXXFLAGS = $(COMP_FLAGS)
build: CPPFLAGS = $(ARMA_FLAGS)
build: LDFLAGS = $(OMP_FLAGS) $(ARMA_LIB)
build: build.o RWGDomain.o fio.o

build.o: CXXFLAGS = $(COMP_FLAGS) $(OMP_FLAGS)
build.o: build.cc

solve: CXXFLAGS = $(COMP_FLAGS)
solve: CPPFLAGS = $(ARMA_FLAGS)
solve: LDLIBS = $(ARMA_LIB)
solve: solve.o RWGDomain.o fio.o

#simpsolve: CXXFLAGS = $(COMP_FLAGS)
#simpsolve: CPPFLAGS = $(ARMA_FLAGS)
#simpsolve: LDLIBS = $(ARMA_LIB)
#simpsolve: simpsolve.o fio.o

integrate.o: point.h integrate.h
RWGDomain.o: RWGDomain.cc point.h RWGDomain.h

clean:
	rm -f *.o *.tri *.mat *.vtk build solve simpsolve

