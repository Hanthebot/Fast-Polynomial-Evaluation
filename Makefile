#*******************************************************************************
#   Project: C++ Library for General Galois Field Arithmetic
#
#   Language: C++ 2007	   
#   Author: Saied H. Khayat
#   Date:   Feb 2013
#   URL: https://github.com/saiedhk
#   
#   Copyright Notice: Free use of this library is permitted under the
#   guidelines and in accordance with the MIT License (MIT).
#   http://opensource.org/licenses/MIT
#
#*******************************************************************************

# GNU Makefile to build the "test_gfelement" executable

CC     = g++
CFLAGS = -std=c++11 -Wall -O3 -o
LFLAGS = -std=c++11 -Wall -O3 -c
GMP = -lgmp -lgmpxx
OBJS   = gfelement.o galoisfield.o modular_arith.o fft.o util.o

fft_finite: fft_finite.cpp fft_finite.h $(OBJS)
	$(CC) $(CFLAGS) fft_finite fft_finite.cpp $(OBJS) $(GMP)

gfelement.o: galois/gfelement.h galois/gfelement.cpp galois/modular_arith.h galois/typedefs.h 
	$(CC) $(LFLAGS) galois/gfelement.cpp $(GMP)

galoisfield.o: galois/galoisfield.h galois/galoisfield.cpp galois/typedefs.h
	$(CC) $(LFLAGS) galois/galoisfield.cpp $(GMP)

modular_arith.o: galois/modular_arith.h galois/modular_arith.cpp galois/typedefs.h
	$(CC) $(LFLAGS) galois/modular_arith.cpp $(GMP)

fft.o: fft.cpp fft_finite.h
	$(CC) $(LFLAGS) fft.cpp $(GMP)

util.o: util.cpp fft_finite.h
	$(CC) $(LFLAGS) util.cpp $(GMP)

clean:
	del *.o

