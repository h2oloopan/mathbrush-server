CXX=g++
CXXFLAGS=-g -I$(LIBMATHRECO_ROOT)/include
LDFLAGS=-L$(LIBMATHRECO_ROOT)/lib

ts : test.o
	g++ -o ts test.o $(LDFLAGS) -static -lmathreco

clean :
	rm -f test.o ts
