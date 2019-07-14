#CPPFLAGS = -I/usr/local/include/eigen3
HOME=/home/crowleyd
CPPFLAGS = -DARMA_DONT_USE_WRAPPER 
CXXFLAGS = -g -O3 -Wall -std=c++11 -fopenmp
LDFLAGS = -fopenmp
LDLIBS = -lstdc++ -lopenblas -lm
CXX=g++
LD=gcc

all: build solve simpsolve

build: build.o RWGDomain.o fio.o
solve: solve.o RWGDomain.o fio.o
simpsolve: simpsolve.o fio.o
proj2: LDLIBS = -lstdc++ -lopenblas -lm -lpthread 
proj2: proj2.o RWGDomain.o integrate.o

spherical.o: spherical.cc spherical.h point.h
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
	rm -f *.o proj *.gnu

