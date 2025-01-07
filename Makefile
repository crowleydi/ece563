
CXX=g++
CC=g++ # for linking

COMP_FLAGS=-O3 -Wall -std=c++14
ARMA_FLAGS=-DARMA_DONT_USE_WRAPPER $(shell pkg-config --cflags armadillo)
ARMA_FLAGS=$(shell pkg-config --cflags armadillo)
BLAS_LIB=$(shell pkg-config --libs armadillo)

ifeq ($(shell uname), Darwin)
OMP_FLAGS=-Xpreprocessor -fopenmp
OMP_LIB=-L/opt/homebrew/opt/libomp/lib -lomp
OMP_INCLUDE=-I/opt/homebrew/opt/libomp/include
else
OMP_FLAGS=-fopenmp
OMP_LIB=-fopenmp
endif

all: build solve simpsolve

build: CXXFLAGS = $(COMP_FLAGS) $(OMP_FLAGS)
build: CPPFLAGS = $(ARMA_FLAGS) $(OMP_INCLUDE)
build: LDFLAGS = $(OMP_LIB) $(BLAS_LIB)
build: build.o RWGDomain.o fio.o

solve: CXXFLAGS = $(COMP_FLAGS)
solve: CPPFLAGS = $(ARMA_FLAGS)
solve: LDLIBS = $(BLAS_LIB)
solve: solve.o RWGDomain.o fio.o

simpsolve: CXXFLAGS = $(COMP_FLAGS)
simpsolve: CPPFLAGS = $(ARMA_FLAGS)
simpsolve: LDLIBS = $(BLAS_LIB)
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

