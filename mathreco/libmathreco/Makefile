CXX=g++
CXXFLAGS=-O3 -I$(LIBMATHRECO_ROOT)/include -DIPAD_RECOGNIZER

objects := $(patsubst %.cc,%.o,$(wildcard *.cc))
objects += $(patsubst %.cpp,%.o,$(wildcard *.cpp))
rangeobjects := $(patsubst %.cpp,%.o,$(wildcard $(LIBMATHRECO_ROOT)/range-expander/*.cpp))
matrixobjects := $(patsubst %.cpp,%.o,$(wildcard $(LIBMATHRECO_ROOT)/matrixreco/*.cpp))

libmathreco.a : $(objects) $(rangeobjects) $(matrixobjects)
	ar -r libmathreco.a $^
	mkdir -p $(LIBMATHRECO_ROOT)/lib
	cp libmathreco.a $(LIBMATHRECO_ROOT)/lib

clean :
	rm -f $(objects) $(rangeobjects) $(matrixobjects) libmathreco.a
