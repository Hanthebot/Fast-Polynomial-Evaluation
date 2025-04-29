#*******************************************************************************
#   Language: C++ 2020
#   
#   Copyright Notice: Free use of this library is permitted under the
#   guidelines and in accordance with the MIT License (MIT).
#   http://opensource.org/licenses/MIT
#
#*******************************************************************************

# GNU Makefile to build the executable

CC     = g++
CFLAGS = -std=c++20 -Wall -O3 -o
LFLAGS = -std=c++20 -Wall -O3 -c
OMP = -fopenmp
GMP = -lgmp -lgmpxx
OBJS   = gfelement.o galoisfield.o modular_arith.o fft.o util.o
OBJS_MULT   = nd_fft.o nd_fft_helper.o fft_brute.o fft_compute.o util_math.o util_multivar.o nd_vector.o

fft_multivar: fft_multivar/fft_multivar.cpp fft_multivar/fft_multivar.h fft_multivar/nd_fft.h fft_multivar/util.h $(OBJS_MULT)
	$(CC) $(OMP) $(CFLAGS) fft_multivar fft_multivar/fft_multivar.cpp $(OBJS_MULT) $(GMP)

nd_fft.o: fft_multivar/nd_fft.cpp fft_multivar/fft_multivar.h fft_multivar/nd_fft.h fft_multivar/util.h
	$(CC) $(OMP) $(LFLAGS) fft_multivar/nd_fft.cpp $(GMP)

nd_fft_helper.o: fft_multivar/nd_fft_helper.cpp fft_multivar/fft_multivar.h fft_multivar/nd_fft.h fft_multivar/util.h
	$(CC) $(LFLAGS) fft_multivar/nd_fft_helper.cpp $(GMP)

fft_brute.o: fft_multivar/fft_brute.cpp fft_multivar/fft_multivar.h
	$(CC) $(LFLAGS) fft_multivar/fft_brute.cpp $(GMP)

fft_compute.o: fft_multivar/fft_compute.cpp fft_multivar/fft_multivar.h fft_multivar/nd_fft.h fft_multivar/util.h
	$(CC) $(OMP) $(LFLAGS) fft_multivar/fft_compute.cpp $(GMP)

util_math.o: fft_multivar/util_math.cpp fft_multivar/fft_multivar.h fft_multivar/util.h
	$(CC) $(LFLAGS) fft_multivar/util_math.cpp $(GMP)

util_multivar.o: fft_multivar/util_multivar.cpp fft_multivar/fft_multivar.h fft_multivar/util.h
	$(CC) $(LFLAGS) fft_multivar/util_multivar.cpp $(GMP)

gfelement.o: galois/gfelement.h galois/gfelement.cpp galois/modular_arith.h galois/typedefs.h 
	$(CC) $(LFLAGS) galois/gfelement.cpp $(GMP)

galoisfield.o: galois/galoisfield.h galois/galoisfield.cpp galois/typedefs.h
	$(CC) $(LFLAGS) galois/galoisfield.cpp $(GMP)

modular_arith.o: galois/modular_arith.h galois/modular_arith.cpp galois/typedefs.h
	$(CC) $(LFLAGS) galois/modular_arith.cpp $(GMP)

nd_vector.o: nd_vector/nd_vector.cpp nd_vector/nd_vector.h
	$(CC) $(OMP) $(LFLAGS) nd_vector/nd_vector.cpp $(GMP)

fft_finite: fft_finite/fft_finite.cpp fft_finite/fft_finite.h $(OBJS)
	$(CC) $(CFLAGS) fft_finite fft_finite/fft_finite.cpp $(OBJS) $(GMP)

fft.o: fft_finite/fft.cpp fft_finite/fft_finite.h
	$(CC) $(LFLAGS) fft_finite/fft.cpp $(GMP)

util.o: fft_finite/util.cpp fft_finite/fft_finite.h
	$(CC) $(LFLAGS) fft_finite/util.cpp $(GMP)

clean:
	del *.o

