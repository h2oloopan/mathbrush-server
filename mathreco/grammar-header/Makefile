CXXFLAGS=-g -I$(LIBMATHRECO_ROOT)/include
LDFLAGS=-L$(LIBMATHRECO_ROOT)/lib

header : header.o
	g++ -o header header.o $(LDFLAGS) -lmathreco

clean :
	rm -f header.o header
