CPPFLAGS= -DNOMATLAB -fopenmp -fPIC -O3

all: c python

phony: c python

# C modules

c: libproxtv.a

libproxtv.a: proxtv.o TVL1opt.o TVL1Wopt.o TVL2opt.o TVLPopt.o TV2Dopt.o TVNDopt.o LPopt.o utils.o condat_fast_tv.o
	ar cr libproxtv.a proxtv.o TVL1opt.o TVL1Wopt.o TVL2opt.o TVLPopt.o TV2Dopt.o TVNDopt.o LPopt.o utils.o condat_fast_tv.o
	ranlib libproxtv.a
	
# Python modules
	
python: proxtv_wrap.cxx _proxtv.so
	
proxtv_wrap.cxx: libproxtv.a proxtv.i proxtv.cc
	swig -c++ -python proxtv.i
	
_proxtv.so: proxtv_wrap.cxx
	python setup.py build_ext --inplace

# Cleaning

clean:
	rm -f *.a *.o *.mexglx *.mexa64 *.so tvlib_wrap.cxx
	rm -rf build
