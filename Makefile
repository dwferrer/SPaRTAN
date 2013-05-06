CXX=icc
CUDACC=nvcc
CXXFLAGS=-O0 -g3 -no-prec-div -openmp -I./ann_float/include/
LIBS := -lgsl -lgslcblas -lm -liomp5 -lcuda -L/usr/local/cuda-5.0/lib64/ -lcudart -L./ann_float/lib/ -lANN

all: testsph libsph.so


testsph:test.cc gravity.o
	icc -fpic -MMD $(CXXFLAGS)  -o$@ $< gravity.o $(LIBS)

libsph.so:pyinterface.o gravity.o
	icc -shared -fpic $(CXXFLAGS) -o $@ $< gravity.o $(LIBS)

pyinterface.o:pyinterface.cc Makefile
	icc -fpic -MMD $(CXXFLAGS) -c $<

gravity.o:gravity.cu Makefile
	nvcc -O3 -Xcompiler -fPIC -c $<

clean:
	rm -f *.o *.d

distclean:
	rm -f *.o *.d *.so testsph

.PHONY: clean distclean
-include pyinterface.d
-include testsph.d
