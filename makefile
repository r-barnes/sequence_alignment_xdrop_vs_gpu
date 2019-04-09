COMPILER = g++
COMPILER_GPU = nvcc
CC = gcc
CFLAGS = -I. -O3 -Wall -Wextra -pedantic -ansi -c -Wno-write-strings
SEQINCLUDE=	-Iseqan

optlist.o:	../optlist/optlist.c ../optlist/optlist.h
	$(CC) $(CFLAGS) $<

test: test-performance.cpp

	$(COMPILER) -std=c++14 $(INCLUDE) -O3 -mavx2 $(SEQINCLUDE) -o test test-performance.cpp

gpu: test-performance-gpu.cu

	$(COMPILER_GPU) -Xcompiler "-std=c++11" $(INCLUDE) $(SEQINCLUDE) test-performance-gpu.cu -o test

gpu_temp: test-performance-gpu-temp.cu

	$(COMPILER_GPU) -std=c++11 $(INCLUDE) $(SEQINCLUDE) test-performance-gpu-temp.cu -o test

clean:
	rm -f *.o
	rm -f test
	rm -f benchmark.txt
	
