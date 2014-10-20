# to see what is inside the library we can use "file mathreco.a" from the command prompt
# use make distclean to clean everything

# general compilation options add -g for debug version
#OPTIONS=-c
#use the following with the app store version
OPTIONS=-c -DNDEBUG

# the device root and sdks
DEVROOT= /Applications/Xcode.app/Contents/Developer/Platforms/iPhoneOS.platform/Developer/
DEVSDKROOT=$(DEVROOT)SDKs/iPhoneOS6.1.sdk
DEVCOMPILER=$(DEVROOT)usr/bin/llvm-gcc-4.2 $(OPTIONS)
DEVAR=$(DEVROOT)usr/bin/ar
DEVRANLIB=$(DEVROOT)usr/bin/ranlib

# the simulator root and sdk and other tools
SIMROOT= /Applications/Xcode.app/Contents/Developer/Platforms/iPhoneSimulator.platform/Developer/
SIMSDKROOT=$(SIMROOT)SDKs/iPhoneSimulator6.1.sdk
SIMCOMPILER=$(SIMROOT)usr/bin/llvm-gcc-4.2 $(OPTIONS)
SIMAR=$(SIMROOT)usr/bin/ar
SIMRANLIB=$(SIMROOT)usr/bin/ranlib

# the MAC compilation
MACCOMPILER=gcc $(OPTIONS)

# the library directories
ROOTDIR=$(CURDIR)
INCLUDEDIR=$(ROOTDIR)/include
LIBDIR=$(ROOTDIR)/lib
SYMBOLSSOURCEDIR=$(ROOTDIR)/libmathreco
MATRIXSOURCEDIR=$(ROOTDIR)/matrixreco
EXPANDERSOURCEDIR=$(ROOTDIR)/range-expander
CASERVERDIR=$(ROOTDIR)/caserver

# the source files
#SYMBOLSSOURCEFILES=parms.cpp annotate.cpp debug.cpp elastic.cpp error.cc feat.cc gridnorm.cc ink-io.cpp memory.cpp norm.cc recog.cc 
#SYMBOLSSOURCEFILES+=segment.cpp grouping.cc verb.cc stroke-alg.cpp symbols.cc memory-policies.cc unix-utils.cpp
#SYMBOLSSOURCEFILES+=MathRecognizer.cpp grammar.cc meas.cpp grammar-values.cc parser.cc bitvec.cc ordered-segments.cc expr-node.cc 
#SYMBOLSSOURCEFILES+=order.cc parser-io.cc relation.cc mathrecognizer-private.cc builder.cc interp.cc match.cc
#MATRIXSOURCEFILES=matrix.cpp grouping-strat.cpp MatrixAnalyzer.cpp matrix-fullyspecified.cpp matrix-underspecified.cpp matrix-graph.cpp
#EXPANDERSOURCEFILES=equiv-transforms.cpp functions.cpp RangeExpander.cpp surrogate-tree.cpp 

SYMBOLSSOURCEFILES=*.cc *.cpp
MATRIXSOURCEFILES=*.cpp
EXPANDERSOURCEFILES=*.cpp
CASERVERFILES=client.cc io.cc log.cc expr.cc CASOperation.cc

# first dependency is the default
TARGETS: SIMTARGET ARM6TARGET ARM7TARGET MACTARGET
	lipo $(LIBDIR)/mathreco_i386.a $(LIBDIR)/mathreco_armv6.a $(LIBDIR)/mathreco_armv7.a -create -output $(LIBDIR)/mathreco.a
	

SIMTARGET:
	cd $(SYMBOLSSOURCEDIR); $(SIMCOMPILER) -I$(INCLUDEDIR) -isysroot $(SIMSDKROOT) -mmacosx-version-min=10.6 -DIPAD_RECOGNIZER -D__IPHONE_OS_VERSION_MIN_REQUIRED=30000 -DDEVVERSION -arch i386 $(SYMBOLSSOURCEFILES)
	cd $(MATRIXSOURCEDIR); $(SIMCOMPILER) -I$(INCLUDEDIR) -isysroot $(SIMSDKROOT) -mmacosx-version-min=10.6 -DIPAD_RECOGNIZER -D__IPHONE_OS_VERSION_MIN_REQUIRED=30000 -DDEVVERSION -arch i386 $(MATRIXSOURCEFILES)
	cd $(EXPANDERSOURCEDIR); $(SIMCOMPILER) -I$(INCLUDEDIR) -isysroot $(SIMSDKROOT) -mmacosx-version-min=10.6 -DIPAD_RECOGNIZER -D__IPHONE_OS_VERSION_MIN_REQUIRED=30000 -DDEVVERSION -arch i386 $(EXPANDERSOURCEFILES)
	cd $(CASERVERDIR); $(SIMCOMPILER) -I$(INCLUDEDIR) -isysroot $(SIMSDKROOT) -mmacosx-version-min=10.6 -DIPAD_RECOGNIZER -D__IPHONE_OS_VERSION_MIN_REQUIRED=30000 -DDEVVERSION -arch i386 $(CASERVERFILES)
	$(SIMAR) rc mathreco_i386.a $(SYMBOLSSOURCEDIR)/*.o $(MATRIXSOURCEDIR)/*.o $(EXPANDERSOURCEDIR)/*.o $(CASERVERDIR)/*.o
	mkdir -p $(LIBDIR)
	$(SIMRANLIB) mathreco_i386.a
	cp mathreco_i386.a $(LIBDIR) 
	make clean

ARM6TARGET:
	cd $(SYMBOLSSOURCEDIR);$(DEVCOMPILER) -I$(INCLUDEDIR) -isysroot $(DEVSDKROOT) -miphoneos-version-min=4.3 -DIPAD_RECOGNIZER -DDEVVERSION -arch armv6 -c $(SYMBOLSSOURCEFILES)
	cd $(MATRIXSOURCEDIR);$(DEVCOMPILER) -I$(INCLUDEDIR) -isysroot $(DEVSDKROOT) -miphoneos-version-min=4.3 -DIPAD_RECOGNIZER -DDEVVERSION -arch armv6 -c $(MATRIXSOURCEFILES)
	cd $(EXPANDERSOURCEDIR);$(DEVCOMPILER) -I$(INCLUDEDIR) -isysroot $(DEVSDKROOT) -miphoneos-version-min=4.3 -DIPAD_RECOGNIZER -DDEVVERSION -arch armv6 -c $(EXPANDERSOURCEFILES)
	cd $(CASERVERDIR);$(DEVCOMPILER) -I$(INCLUDEDIR) -isysroot $(DEVSDKROOT) -miphoneos-version-min=4.3 -DIPAD_RECOGNIZER -DDEVVERSION -arch armv6 -c $(CASERVERFILES)
	$(DEVAR) rc mathreco_armv6.a $(SYMBOLSSOURCEDIR)/*.o $(MATRIXSOURCEDIR)/*.o $(EXPANDERSOURCEDIR)/*.o $(CASERVERDIR)/*.o
	$(DEVRANLIB) mathreco_armv6.a
	cp mathreco_armv6.a $(LIBDIR)
	make clean

ARM7TARGET:	
	cd $(SYMBOLSSOURCEDIR);$(DEVCOMPILER) -I$(INCLUDEDIR) -isysroot $(DEVSDKROOT) -miphoneos-version-min=4.3 -DIPAD_RECOGNIZER -DDEVVERSION -arch armv7 -c $(SYMBOLSSOURCEFILES)
	cd $(MATRIXSOURCEDIR);$(DEVCOMPILER) -I$(INCLUDEDIR) -isysroot $(DEVSDKROOT) -miphoneos-version-min=4.3 -DIPAD_RECOGNIZER -DDEVVERSION -arch armv7 -c $(MATRIXSOURCEFILES)
	cd $(EXPANDERSOURCEDIR);$(DEVCOMPILER) -I$(INCLUDEDIR) -isysroot $(DEVSDKROOT) -miphoneos-version-min=4.3 -DIPAD_RECOGNIZER -DDEVVERSION -arch armv7 -c $(EXPANDERSOURCEFILES)
	cd $(CASERVERDIR);$(DEVCOMPILER) -I$(INCLUDEDIR) -isysroot $(DEVSDKROOT) -miphoneos-version-min=4.3 -DIPAD_RECOGNIZER -DDEVVERSION -arch armv7 -c $(CASERVERFILES)
	$(DEVAR) rc mathreco_armv7.a $(SYMBOLSSOURCEDIR)/*.o $(MATRIXSOURCEDIR)/*.o $(EXPANDERSOURCEDIR)/*.o $(CASERVERDIR)/*.o
	$(DEVRANLIB) mathreco_armv7.a
	cp mathreco_armv7.a $(LIBDIR)
	make clean

MACTARGET:
	cd $(SYMBOLSSOURCEDIR); $(MACCOMPILER) -I$(INCLUDEDIR) -DDEVVERSION -arch i386 -c $(SYMBOLSSOURCEFILES)
	cd $(MATRIXSOURCEDIR); $(MACCOMPILER) -I$(INCLUDEDIR) -DDEVVERSION -arch i386 -c $(MATRIXSOURCEFILES)
	cd $(EXPANDERSOURCEDIR); $(MACCOMPILER) -I$(INCLUDEDIR) -DDEVVERSION -arch i386 -c $(EXPANDERSOURCEFILES)
	cd $(CASERVERDIR); $(MACCOMPILER) -I$(INCLUDEDIR) -DDEVVERSION -arch i386 -c $(CASERVERFILES)
	ar -r libmathreco.a $(SYMBOLSSOURCEDIR)/*.o $(MATRIXSOURCEDIR)/*.o $(EXPANDERSOURCEDIR)/*.o $(CASERVERDIR)/*.o
	cp libmathreco.a $(LIBDIR) 
	make clean
	
ARMTARGETS:ARM6TARGET ARM7TARGET
	lipo $(LIBDIR)/mathreco_armv6.a $(LIBDIR)/mathreco_armv7.a -create -output $(LIBDIR)/mathrecoARM.a

clean:
	rm -f $(SYMBOLSSOURCEDIR)/*.o $(MATRIXSOURCEDIR)/*.o $(EXPANDERSOURCEDIR)/*.o $(CASERVERDIR)/*.o *.a
	
distclean:
	make clean
	rm -f $(LIBDIR)/*.*
