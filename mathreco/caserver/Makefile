CFLAGS=-g -I$(LIBMATHRECO_ROOT)/include
CXXFLAGS=-g -I$(LIBMATHRECO_ROOT)/include

PYTHON=python2.7

caserver: server.o log.o expr.o io.o
	g++ -o caserver $^ -L$(LIBMATHRECO_ROOT)/lib -lmathreco -I/usr/include/$(PYTHON) -L/usr/lib/$(PYTHON)/config -Xlinker -l$(PYTHON) -lm -lutil

server.o: server.cc
	g++ $(CXXFLAGS) -I/usr/include/$(PYTHON) -c -o server.o $^ 

clientobjs := client.o io.o expr.o log.o CASOperation.o

libclient: $(clientobjs)
	ar -r libclient.a $^
	cp libclient.a $(LIBMATHRECO_ROOT)/lib

caclient: libclient iclient.o
	g++ -o caclient iclient.o -L$(LIBMATHRECO_ROOT)/lib -lmathreco -lclient

clean:
	rm -f caserver caclient libclient.a *.o
