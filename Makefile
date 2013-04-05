CXX=icc
CUDACC=nvcc
CXXFLAGS=-O3 -openmp 
LIBS := -lgsl -lgslcblas -lm -liomp5 -lcuda -L/usr/local/cuda-5.0/lib64/ -lcudart

all: testsph libsph.so


testsph:test.cc gravity.o
	icc -fpic $(CXXFLAGS) -o$@ $^ $(LIBS)

libsph.so:gravity.o sph.o
	icc -shared -fpic $(CXXFLAGS) -o $@ $< $(LIBS)

sph.o:sph.cc
	icc -fpic -MMD $(CXXFLAGS) -c $<

gravity.o:gravity.cu
	nvcc --use_fast_math -Xcompiler -fPIC -c $<


clean:
	rm -f *.o *.d

distclean:
	rm -f *.o *.d *.so testsph

.PHONY: clean distclean
-include sph.d
