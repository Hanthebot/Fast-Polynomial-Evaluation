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
GMP = -lgmp -lgmpxx
OBJS   = gfelement.o galoisfield.o modular_arith.o fft.o util.o
OBJS_MULT   = gfelement.o galoisfield.o modular_arith.o nd_fft.o util_multivar.o nd_vector.o

fft_finite: fft_finite/fft_finite.cpp fft_finite/fft_finite.h $(OBJS)
	$(CC) $(CFLAGS) fft_finite fft_finite/fft_finite.cpp $(OBJS) $(GMP)

fft_multivar: fft_multivar/fft_multivar.cpp fft_multivar/fft_multivar.h $(OBJS_MULT)
	$(CC) $(CFLAGS) fft_multivar fft_multivar/fft_multivar.cpp $(OBJS_MULT) $(GMP)

gfelement.o: galois/gfelement.h galois/gfelement.cpp galois/modular_arith.h galois/typedefs.h 
	$(CC) $(LFLAGS) galois/gfelement.cpp $(GMP)

galoisfield.o: galois/galoisfield.h galois/galoisfield.cpp galois/typedefs.h
	$(CC) $(LFLAGS) galois/galoisfield.cpp $(GMP)

modular_arith.o: galois/modular_arith.h galois/modular_arith.cpp galois/typedefs.h
	$(CC) $(LFLAGS) galois/modular_arith.cpp $(GMP)

fft.o: fft_finite/fft.cpp fft_finite/fft_finite.h
	$(CC) $(LFLAGS) fft_finite/fft.cpp $(GMP)

nd_fft.o: fft_multivar/nd_fft.cpp fft_multivar/fft_multivar.h
	$(CC) $(LFLAGS) fft_multivar/nd_fft.cpp $(GMP)

util.o: fft_finite/util.cpp fft_finite/fft_finite.h
	$(CC) $(LFLAGS) fft_finite/util.cpp $(GMP)

util_multivar.o: fft_multivar/util_multivar.cpp fft_multivar/fft_multivar.h
	$(CC) $(LFLAGS) fft_multivar/util_multivar.cpp $(GMP)

nd_vector.o: nd_vector/nd_vector.cpp nd_vector/nd_vector.h
	$(CC) $(LFLAGS) nd_vector/nd_vector.cpp $(GMP)

clean:
	del *.o

